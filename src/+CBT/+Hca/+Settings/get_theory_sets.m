function [ sets ] = get_theory_sets( sets,n1)
 


prompt = cell(n1+1, 1);
     prompt{1} = sprintf('Enter conc. of DNA :', 0.2);
     definput{1}='0.2'
for i = 2:2
    prompt{i} = sprintf('Enter conc. of ligand %d:', i-1);
     definput{i}='6'
end

for i = 3:n1+1
    prompt{i} = sprintf('Enter conc. of ligand %d:', i-1);
     definput{i}='0.02'
end


dlgTitle = 'Input Dialog';
dims = [1 50];  
answer = inputdlg(prompt, dlgTitle, dims,definput);
        for i=2:n1+1
       % ligandconc(i)  = str2double(answer{i});
        sets.ligandconc(i-1)= str2double(answer{i});
        end 
 
        if ~isempty(answer)
        sets.concDNA= str2double(answer{1});
        sets.concN=str2double(answer{2});
        sets.concY=str2double(answer{3});
        else
               disp('Default theory settings are being used');
         end
 
        
