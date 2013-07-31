%Author: Rahul Mohan, 2013

%Once you align your reads using Bowtie,TopHat you will get an output in the .sam
%format. Convert it to a txt (duplicate it and just change the extension
%from .sam to .txt) and load it into Excel. The fifth column contains map
%scores of each base-call:
%
%255 = unique mapping
%3 = maps to 2 locations in the genome
%2 = maps to 3 locations in the genome
%1 = maps to 4-9 locations in the genome
%0 = maps to 10 locations in the genome
%
%Extract that column; this program calculates how accurate Bowtie was in
%aligning to the reference genome with that as an input


%%%% Input your vector in variable MapScores %%%%

MapScores = xlsread('TopHatBrainOutput.xlsx'); %Change 'MapScoresBowtie.xlsx' to whatever file you saved containing vector of map scores

%Define error rates of each possible map score

P_3 = 1/2; %Mapped to 2 locations in genome, so error rate is 1/2
P_2 = 2/3; %Mapped to 3 locations in genome, so error rate is 2/3
P_1 = mean([3/4,4/5,5/6,7/8,8/9]); %Mapped to 4-9 locations in genome, so error rate is %average of [3/4,4/5,5/6,7/8,8/9]
P_0 = 9/10; %Mapped to 10 locations in genome, so error rate is 9/10

%255 is not included because it means that the alignment was perfectly
%accurate

%Phred calculations

phred255 = 60; %Basically a Phred of 60 is considered to be the lowest probability of error there is, so assign it to 255 (unique mapping)
phred3 = -10 * log10(P_3);
phred2 = -10 * log10(P_2);
phred1 = -10 * log10(P_1);
phred0 = -10 * log10(P_0);

nonuniquePairs = 0;
error = zeros(size(MapScores));
totalSequences = length(MapScores);

%Replace 255,3,2,1,0 with proper phred score calculations: Preprocessing step for loop

MapScores(MapScores==255) = phred255;
MapScores(MapScores==3) = phred3;
MapScores(MapScores==2) = phred2;  %Precalculated Phred Scores
MapScores(MapScores==1) = phred1;
MapScores(MapScores==0) = phred0;


%Loop to count the number of different map scores and phred error


for i = 1:length(MapScores)
    
    error(i) = double((10 .^ ((-MapScores(i))/10)) ); %Phred Prob. of Error Formula
    
    if MapScores(i) ~= 60
        nonuniquePairs = nonuniquePairs + 1;
    end
    
end

%%% What the above loop says is that: add the corresponding phred_error probability to a
%%% vector of zeros at that specific index.

%Specific Phred Errors

phred3_error = ((10 .^ ((-phred3)/10)) ); 
phred2_error = ((10 .^ ((-phred2)/10)) );
phred1_error = ((10 .^ ((-phred1)/10)) );
phred0_error = ((10 .^ ((-phred0)/10)) );

%Phred Score Calculations
         
probError = mean(error); %Average of the error vector

%% printf's which calculate final output information
         
fprintf(1,'\nNumber of Map Scores of 255: %d \n\n',length(MapScores(MapScores==60)));

fprintf(1,'\nNumber of Map Scores of 3: %d \n\n',length(error(error==phred3_error)));

fprintf(1,'\nNumber of Map Scores of 2: %d \n\n',length(error(error==phred2_error)));

fprintf(1,'\nNumber of Map Scores of 1: %d \n\n',length(error(error==phred1_error)));

fprintf(1,'\nNumber of Map Scores of 0: %d \n\n',length(error(error==phred0_error)));

fprintf(1,'\nNumber of non-unique mappings: %d \n\n',nonuniquePairs);

fprintf(1,'\n%d sequences with an average alignment error of: ~ %d \n\n',totalSequences,(probError));
