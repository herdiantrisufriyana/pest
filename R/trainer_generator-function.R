#' Create classification training function for DI-VNN given lambda
#'
#' Create training function for deep-insight visible neural network given 
#' lambda (L2-norm regularization parameter) for classification task.
#'
#' @param tidy_set Tidy set, an ExpressionSet that includes training set.
#' @param path Path and the prefix, a character of the path that includes the 
#' prefix for any result saved from this training.
#' @param class_weight Outcome weight, a vector of numeric for nonevent and 
#' event.
#' @param epochs Maximum epoch, an integer of the maximum epoch.
#' @param patience Minimum early stop, an integer of the minimum epoch that 
#' early stop is applied.
#' @param batch_size Batch size, an integer of the number of instances used 
#' for each step or iteration.
#' @param warm_up Warm-up steps, a numeric of proportion of the early training 
#' instances during warming up learning rate from the one thirty-second of it.
#' @param lr Initial learning rate, a numeric of the initial learning rate.
#' @param min_lr Minimum learning rate, a numeric of the minimum learning rate 
#' can be applied after reduced when the evaluation are getting plateau.
#' @param tuning_mode Tuning mode, a logical whether the returning output is 
#' one for hyperparameter tuning. See Return.
#' @param checkpoint Checkpoint, a logical whether the checkpoint of the model 
#' used to continue the training.
#' @param verbose Verbosity, an integer of 0 or 1, that is respectively false 
#' or true whether the training messages are shown.
#'
#' @details
#' This function is used to create another function. It will need lambda value 
#' as the only argument. Lambda is L2-norm regularization parameter.
#' 
#' @return The output depends on \code{tuning_mode}. If true, the output is the 
#' evaluation score, a numeric of the score. If false, the output is a list of 
#' Keras model object, early-stopping callback object, training index vector, 
#' validation index vector, a list of history of the Keras training, and a list  
#' of resampled evaluation results.
#' 

trainer_generator=function(tidy_set
                           ,path=NULL
                           ,class_weight=NULL
                           ,epochs=500
                           ,patience=round(500/2)
                           ,batch_size=128
                           ,warm_up=0.05
                           ,lr=2^-6
                           ,min_lr=2^-6/32
                           ,tuning_mode=F
                           ,checkpoint=F
                           ,verbose=1){
  
  # This function builds another function.
  function(lambda){
    
    # Minimize memory accumulation.
    refresh_session()
    
    # Construct the DI-VNN architecture of the Keras model.
    ontonet=
      tidy_set %>%
      generator.ontonet(
        path=NULL
        ,l2_norm=lambda
        ,output_unit=1
        ,output_activation='sigmoid'
      )
    
    # If checkpoint is true, then check existing checkpoint in the path.
    # If existed, then load the model weights 
    # from the previous (usually interrupted) training.
    if(checkpoint){
      if(dir.exists(path)){
        cat('Checkpoint is detected. Restore checkpoint: ',sep='')
        cat('...',sep='')
        load_model_weights_tf(ontonet,paste0(path,'/ckpt'))
        cat(' Done.\n',sep='')
      }
    }
    
    # Compile the Keras model.
    ontonet %>%
      compile(
        optimizer=optimizer_sgd(lr=lr,momentum=0.9)
        ,loss='mean_squared_error'
        ,loss_weights=
          c(rep(0.3/(0.3*(length(.$outputs)-1)+1),length(.$outputs)-1)
            ,1/(0.3*(length(.$outputs)-1)+1))
        ,metrics=c(
          tf$keras$metrics$AUC(name='roc')
          ,tf$keras$metrics$TruePositives(name='tp')
          ,tf$keras$metrics$FalseNegatives(name='fn')
          ,tf$keras$metrics$FalsePositives(name='fp')
          ,tf$keras$metrics$TrueNegatives(name='tn')
        )
      )
    
    # Set early-stopping criteria.
    cb_early_stopping=
      callback_early_stopping(
        monitor='val_root_roc'
        ,mode='max'
        ,min_delta=0.001
        ,patience=patience
        ,restore_best_weights=T
      )
    
    # Set learning-rate reduction criteria.
    cb_lr_reduction=
      callback_reduce_lr_on_plateau(
        monitor='val_root_roc'
        ,factor=0.96
        ,patience=1
        ,verbose=0
        ,mode='max'
        ,min_delta=0.01
        ,cooldown=0
        ,min_lr=min_lr
      )
    
    # Build a function used at the beginning of epoch.
    on_epoch_begin_fun=function(epoch,logs){
      current_epoch_lr=k_get_value(ontonet$optimizer$lr)
      lr_gradient<<-seq(current_epoch_lr/32,current_epoch_lr,len=warm_up_steps)
    }
    
    # Build a function used at the end of training step or iteration.
    on_train_batch_end_fun=function(batch,logs){
      if((batch+1)<=warm_up_steps){
        k_set_value(ontonet$optimizer$lr,lr_gradient[batch+1])
      }
      cat(sprintf('\r'),sep='')
      cat('Step: ',batch+1,'/',steps_per_epoch,'; ',sep='')
      cat('LR: ',round(k_get_value(ontonet$optimizer$lr),10),'; ',sep='')
      cat('Total: ',append=T,sep='')
      cat('loss=',round(logs[["loss"]],4),'; ',sep='')
      cat('Root: ',append=T,sep='')
      cat('loss=',round(logs[["root_loss"]],4),'; ',sep='')
      cat('ROC=',round(logs[["root_roc"]],4),'; ',sep='')
      cat('TP=',round(logs[["root_tp"]]),'; ',sep='')
      cat('FN=',round(logs[["root_fn"]]),'; ',sep='')
      cat('FP=',round(logs[["root_fp"]]),'; ',sep='')
      cat('TN=',round(logs[["root_tn"]]),'     ',sep='')
    }
    
    # Build a function used at the end of epoch 
    # after the in-epoch validation.
    on_epoch_end_fun=function(epoch,logs){
      cat('\n',sep='')
      cat('Epoch: ',epoch+1,'/',epochs,'; ',sep='')
      cat('LR: ',round(k_get_value(ontonet$optimizer$lr),10),'; ',sep='')
      cat('Total: ',sep='')
      cat('loss=',round(logs[["val_loss"]],4),'; ',sep='')
      cat('Root: ',append=T,sep='')
      cat('loss=',round(logs[["val_root_loss"]],4),'; ',sep='')
      cat('ROC=',round(logs[["val_root_roc"]],4),'; ',sep='')
      cat('TP=',round(logs[["val_root_tp"]]),'; ',sep='')
      cat('FN=',round(logs[["val_root_fn"]]),'; ',sep='')
      cat('FP=',round(logs[["val_root_fp"]]),'; ',sep='')
      cat('TN=',round(logs[["val_root_tn"]]),sep='')
      cat('\n')
    }
    
    # Build a function used at the end of validation step
    # for each time of resampling.
    on_test_batch_end_fun=function(batch,logs){
      cat(sprintf('\r'),sep='')
      cat('Evaluate: ',batch+1,'/',steps,'; ',sep='')
      cat('LR: ',round(k_get_value(ontonet$optimizer$lr),10),'; ',sep='')
      cat('Total: ',append=T,sep='')
      cat('loss=',round(logs[["loss"]],4),'; ',sep='')
      cat('Root: ',append=T,sep='')
      cat('loss=',round(logs[["root_loss"]],4),'; ',sep='')
      cat('ROC=',round(logs[["root_roc"]],4),'; ',sep='')
      cat('TP=',round(logs[["root_tp"]]),'; ',sep='')
      cat('FN=',round(logs[["root_fn"]]),'; ',sep='')
      cat('FP=',round(logs[["root_fp"]]),'; ',sep='')
      cat('TN=',round(logs[["root_tn"]]),'     ',sep='')
    }
    
    # Build a function used at the end of validation.
    on_test_end_fun=function(logs){
      cat('\n')
    }
    
    # If verbosity is true, then use the functions at the multiple offsets. 
    # Otherwise, only use early-stopping and learning-rate reduction functions.
    if(verbose>0){
      callback_list=
        c(cb_early_stopping
          ,cb_lr_reduction
          ,callback_lambda(
            on_epoch_begin=on_epoch_begin_fun
            ,on_train_batch_end=on_train_batch_end_fun
            ,on_epoch_end=on_epoch_end_fun
          ))
    }else{
      callback_list=c(cb_early_stopping,cb_lr_reduction)
    }
    
    # if checkpoint is true, then save a checkpoint at the end of each epoch.
    if(checkpoint){
      cb_checkpoint=
        callback_model_checkpoint(
          paste0(path,'/ckpt')
          ,monitor='val_root_loss'
          ,verbose=0
          ,save_best_only=T
          ,save_weights_only=T
          ,mode='min'
          ,save_freq='epoch'
        )
      callback_list=c(callback_list,cb_checkpoint)
    }
    
    # Split data for in-epoch validation.
    suppressWarnings(set.seed(33,sample.kind=sample.kind))
    index=sample(1:dim(tidy_set)[2],dim(tidy_set)[2],F)
    
    val_i=
      index %>%
      sample(round(0.2*length(.)),F)
    
    train_i=index[!index%in%val_i]
    
    # Determine steps per epoch for training and validation, 
    # including the number of steps for warming up
    steps_per_epoch=ceiling(length(train_i)/batch_size)
    validation_steps=ceiling(length(val_i)/batch_size)
    warm_up_steps=ceiling(steps_per_epoch*warm_up)
    
    # Assign outcome weights as a vector as much as the training instances.
    sample_weights=
      tidy_set %>%
      phenoData() %>%
      pData() %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        data.frame(sample_weight=class_weight) %>%
          rownames_to_column(var='outcome') %>%
          mutate(outcome=as.integer(outcome))
        ,by='outcome'
      ) %>%
      lapply(X=1,Y=.,function(X,Y){
        setNames(Y$sample_weight,Y$id)
      }) %>%
      .[[1]]
    
    # Train the model.
    cat('Warm-up: ',warm_up_steps,'steps\n')
    history=
      ontonet %>%
      fit_generator(
        generator=
          generator.ontoarray(
            tidy_set
            ,train_i
            ,batch_size=batch_size
            ,sample_weights=sample_weights
          )
        ,steps_per_epoch=steps_per_epoch
        ,validation_data=
          generator.ontoarray(
            tidy_set
            ,val_i
            ,batch_size=batch_size
            ,sample_weights=sample_weights
          )
        ,validation_steps=validation_steps
        ,class_weight=
          ontonet %>%
          get_config() %>%
          .$layers %>%
          sapply(function(x)x$name) %>%
          .[!str_detect(.,'_')] %>%
          .[.!='ontoarray'] %>%
          lapply(X=1,Y=.,Z=class_weight,function(X,Y,Z){
            lapply(X=seq(length(Y)),Y=Z,function(X,Y){
              Y
            }) %>%
              setNames(Y)
          }) %>%
          .[[1]]
        ,epochs=epochs
        ,callbacks=callback_list
        ,view_metrics=F
        ,verbose=0
      )
    
    # Prepare validation for number of steps, bootstraps, and the callbacks.
    suppressWarnings(set.seed(33,sample.kind=sample.kind))
    steps=ceiling(dim(tidy_set)[2]/dim(tidy_set)[2])
    steps=ceiling(length(val_i)/batch_size)
    boots=30
    if(verbose>0){
      if(tuning_mode){
        eval_callback_list=c(cb_early_stopping,cb_lr_reduction)
      }else{
        eval_callback_list=
          c(cb_early_stopping
            ,cb_lr_reduction
            ,callback_lambda(
              on_test_batch_end=on_test_batch_end_fun
              ,on_test_end=on_test_end_fun
            ))
      }
    }else{
      eval_callback_list=c(cb_early_stopping,cb_lr_reduction)
    }
    
    # Evaluate the model using the validation split.
    evaluation=
      ontonet %>%
      lapply(X=1:boots,Y=.,Z=tidy_set,function(X,Y,Z){
        if(verbose>0){
          cat(sprintf('\r'),'Evaluation ',X,'/',boots,' ',sep='')
        }
        index=sample(val_i,length(val_i),T)
        ontonet %>%
          evaluate_generator(
            generator=generator.ontoarray(
              Z
              ,index
              ,batch_size=batch_size
              ,sample_weights=sample_weights
            )
            ,steps=steps
            ,callbacks=eval_callback_list
          )
      })
    
    # Return output based on tuning mode setting.
    if(tuning_mode){
      score=
        evaluation %>%
        sapply(function(x)x$root_roc) %>%
        mean()
      cat('ROC=',round(score,4),'\n',sep='')
      return(list(Score=score))
    }else{
      list(
        ontonet=ontonet
        ,early_stopping=cb_early_stopping
        ,train_index=train_i
        ,validation_index=val_i
        ,history=history
        ,evaluation=evaluation
      )
    }
  }
}