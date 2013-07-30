%Once you align your reads using BWA, you will get an output in the .sam
%format. Convert it to a txt (duplicate it and just change the extension
%from .sam to .txt) and load it into Excel. The fifth column contains phred
%quality scores of each sequence(row) in the .sam file.
%
%As the phred quality score goes up the probability of an error at that
%sequence goes down. For example, a phred score of 60 means that there is a
%1/1000000 chance that the sequence is aligned incorrectly.
%
%Extract that column; this program calculates how accurate your overall
%aligned reads are from BWA.


%%%% Input your vector in MapScores %%%%

MapScores = xlsread('MapScoresBWA.xlsx'); %Change 'MapScoresBowtie.xlsx' to whatever file you saved containing vector of map scores

%Initialize variables

threshold = 30; % Standard phred threshold, scores below 30 are considered bad quality %reads, vice versa; change if you want
bad_reads = 0;
good_reads = 0;
totalSequences = length(MapScores);

error = zeros(size(MapScores)); %Initialize error vector

%Preprocessing step is not needed like the TopHat and Bowtie calculator because
%the input vector here, outputted by BWA,already contains the precise phred values

%Loop to calculate bad reads, good reads, and alignment error

for i = 1:length(MapScores)
    
    error(i) = double((10 .^ ((-MapScores(i))/10)) ); %Phred Prob. of Error Formula
    
    if MapScores(i) < threshold
        bad_reads = bad_reads + 1;
    end
    
    if MapScores(i) >= threshold
        good_reads = good_reads + 1;
    end
      
end

%Phred Score Calculations

probError = mean(error); %Average of error vector
         
fprintf(1,'\nNumber of Bad Map Scores: %d \n\n', bad_reads);

fprintf(1,'\nNumber of Good Map Scores: %d \n\n',good_reads);

fprintf(1,'\n%d sequences with an average alignment error of: ~ %d \n\n',totalSequences,(probError));
