function [filesBpsPerPixel, errorMsg] = prompt_files_bps_per_pixel(filepaths, hPanel)
    hButtonContinue = [];
    filesBpsPerPixel = [];
    errorMsg = [];
    if isempty(filepaths)
        return;
    end

    function onCompleteCallback(outputFilesBpsPerPixel)
        filesBpsPerPixel = outputFilesBpsPerPixel;
        delete(hButtonContinue);
    end
    import OldDBM.General.Import.async_prompt_files_bps_per_pixel;
    hButtonContinue = async_prompt_files_bps_per_pixel(filepaths, hPanel, @onCompleteCallback);
    waitfor(hButtonContinue);
    if isempty(filesBpsPerPixel) || any(isnan(filesBpsPerPixel))
        errorMsg = 'Valid inputs were not provided for bps/pixel';
    end
end