%% Define path and open EEGlab
clc
clear all
close all
lm_Conf                 = struct();
lm_Conf.path            = '~/Documentos/Repos/CuBaPeTo/';
lm_Conf.customFunsPath  = ['cstFuns/'];     
lm_Conf.csvPath         = ['csv_new/']; 

addpath(genpath(lm_Conf.path))
addpath(genpath(lm_Conf.customFunsPath))

% [lm_Conf, SUJ]= defino_path(lm_Conf);

lm_Conf.rFunctionsPath        = [lm_Conf.path  '/R_functions/'];
lm_Conf.bashPath              = [lm_Conf.path  '/bash_functions/'];

%% Run LMM with 500 permutations across both ranEf

fixEf   = 'freq + pred:tipo + palnum:tipo';
ranEf   = '(1|pal) + (1|suj_id)';
nIter   = 500;
modType = 'lmm';
nCores  = 4;
lm_Conf.nTimes                = 103;
lm_Conf.lmmOutPath            = ['lmm_results_CuBaPeTo/results/'];
lm_Conf.nohupOutPath          = ['lmm_results_CuBaPeTo/nohupOutputs/'];
lm_Conf.permutationMatPath    = ['lmm_results_CuBaPeTo/permutations/']; 
lm_Conf.customFunsPath        = ['cstFuns/'];     

lm_Conf.permutationVariable   = 'across';

lm_parallelRunLMM(fixEf, ranEf, nIter, modType, nCores, lm_Conf)
