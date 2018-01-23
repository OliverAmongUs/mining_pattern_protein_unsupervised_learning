%% 5 sample test
clear
proteinARGs = cell(0);
proteinARGs{end+1} = GenerateProteinARG('protein/test1/1ojyA03.csv');
proteinARGs{end+1} = GenerateProteinARG('protein/test1/1ojyB03.csv');
proteinARGs{end+1} = GenerateProteinARG('protein/test4/1rlhA01.csv');
proteinARGs{end+1} = GenerateProteinARG('protein/test5/4glnE00.csv');
%proteinARGs{end+1} = GenerateProteinARG('protein/test6/1hknD00.csv');
mdl = sprMDL(proteinARGs, 2);
original_result = zeros([1,length(proteinARGs)]);
original_score = zeros([1,length(proteinARGs)]);
for i = 1:length(proteinARGs)
    [result, score] = mdl.checkPattern(proteinARGs{i});
    original_result(i) = result;
    original_score(i) = score;
end
original_detect_rate = sum(original_result)/length(original_score);