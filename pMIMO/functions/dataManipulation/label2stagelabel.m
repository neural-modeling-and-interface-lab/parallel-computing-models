function [stageLabel, wildcard] = label2stagelabel(label,stage)
    switch stage
        case 1
            wildcard = 'Path';
        case 2 
            wildcard = 'ReMLE';
        case 3
            wildcard = 'CVCrit';
    end
 stageLabel = [label wildcard];