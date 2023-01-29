%% Code Explaination must be provided. 
function threedfile = twodvect3dmat(twodfile_array)
    N = size(twodfile_array,1); %Columns in the twodfile 
    M = size(twodfile_array,2); %Rows in the twodfile
    %M_1 = 10^(log10(M)/2); 
    M_1 = sqrt(M);
    twodfile_array_trans = transpose(twodfile_array);
    threedfile = pagetranspose(reshape(twodfile_array_trans,M_1,M_1,N));
end