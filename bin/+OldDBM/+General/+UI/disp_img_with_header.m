function [] = disp_img_with_header(hAxis, imgArr, imgHeaderText)
    imagesc(imgArr, 'Parent', hAxis);
    hold(hAxis, 'on');
    colormap(hAxis, gray());
    
    import OldDBM.General.UI.set_centered_header_text;
    set_centered_header_text(hAxis, imgHeaderText);
    hold(hAxis, 'off');
end