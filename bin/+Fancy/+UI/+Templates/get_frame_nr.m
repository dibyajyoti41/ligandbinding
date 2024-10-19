function [ timeFramesNr] = get_frame_nr(promptForTimeFr,default_timeFramesNr,titleText)
   % get number of timeframes from user
	timeFramesNr = default_timeFramesNr;
    
    if promptForTimeFr ~= 0
        options.Resize = 'on';
        % sets number of time-frames
        answer = inputdlg(...
            { ...
                'Choose number of time-frames', .....
            }, ... % prompt
            titleText, ... % dialog title
            1, ... % number of lines
            { ...
                num2str(timeFramesNr), ...
            },options);
        if ~isempty(answer)
            timeFramesNr = str2double(answer{1});
        end

        isAcceptable = struct;
        isAcceptable.timeFramesNr = (timeFramesNr >= 0);

        if not(isAcceptable.timeFramesNr)
            timeFramesNr = default_timeFramesNr;
            warning('Bad input for nr of timeframes! Try again!');
        end        
    end
end