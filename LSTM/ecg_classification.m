%% Classify ECG Signals Using Long Short-Term Memory Networks
% Based on https://de.mathworks.com/help/signal/ug/classify-ecg-signals-using-long-short-term-memory-networks.html

%% Setup
addpath('/Users/nguyentiendung/Documents/MATLAB/Examples/R2022b/deeplearning_shared/ClassifyECGSignalsUsingLSTMNetworksWithGPUAccelerationExample')

%% Load data
if ~isfile('PhysionetData.mat')
    ReadPhysionetData         
end
load PhysionetData

%% Examine data
L = cellfun(@length,Signals);
h = histogram(L);
xticks(0:3000:18000);
xticklabels(0:3000:18000);
title('Signal Lengths')
xlabel('Length')
ylabel('Count')

normal = Signals{1};
aFib = Signals{4};

subplot(2,1,1)
plot(normal)
title('Normal Rhythm')
xlim([4000,5200])
ylabel('Amplitude (mV)')
text(4330,150,'P','HorizontalAlignment','center')
text(4370,850,'QRS','HorizontalAlignment','center')

subplot(2,1,2)
plot(aFib)
title('Atrial Fibrillation')
xlim([4000,5200])
xlabel('Samples')
ylabel('Amplitude (mV)')

%% Prepare data for training
% make all signals in the input array 9000 samples long
[Signals,Labels] = segmentSignals(Signals,Labels); 

%% Oversampling because the ratio of AF to normal signals is 1:7
afibX = Signals(Labels == 'A');
afibY = Labels(Labels == 'A');

normalX = Signals(Labels == 'N');
normalY = Labels(Labels == 'N');

% train-test split
[trainIndA,~,testIndA] = dividerand(size(afibX, 1),0.9,0.0,0.1);
[trainIndN,~,testIndN] = dividerand(size(normalX, 1),0.9,0.0,0.1);

XTrainA = afibX(trainIndA);
YTrainA = afibY(trainIndA);

XTrainN = normalX(trainIndN);
YTrainN = normalY(trainIndN);

XTestA = afibX(testIndA);
YTestA = afibY(testIndA);

XTestN = normalX(testIndN);
YTestN = normalY(testIndN);

% oversampling
XTrain = [repmat(XTrainA(1:634),7,1); XTrainN(1:4438)];
YTrain = [repmat(YTrainA(1:634),7,1); YTrainN(1:4438)];

XTest = [repmat(XTestA(1:70),7,1); XTestN(1:490)];
YTest = [repmat(YTestA(1:70),7,1); YTestN(1:490);];

%% Define LSTM network architecture
layers = [ ...
    sequenceInputLayer(1)
    bilstmLayer(100, 'OutputMode','last')
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer
    ]

options = trainingOptions('adam', ...
    'MaxEpochs',10, ...
    'MiniBatchSize', 150, ...
    'InitialLearnRate', 0.01, ...
    'SequenceLength', 1000, ...
    'GradientThreshold', 1, ...
    'ExecutionEnvironment',"auto",...
    'plots','training-progress', ...
    'Verbose',false);

%% Train the network
ecgnet = trainNetwork(XTrain,YTrain,layers,options);
cd "C:\Users\tien-\OneDrive\Documents\Studium\Research\BrainData\matlab\Playground\LSTM"
save ecgnet

%% Load network after training
load ecgnet

%% Classification
trainPred = classify(ecgnet,XTrain,"SequenceLength",1000);