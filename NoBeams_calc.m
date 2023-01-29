% Number of beams required to Compute the RCS upto 1% error level 
function [N_arr,beam_array_eps1_mnst_arr,err] = NoBeams_calc(Z_r)
eps1 = logspace(-10,0,1000);
a = 1;
[row_elem,colmn_elem] = size(Z_r);
beam_array_eps1_mnst_arr = zeros(colmn_elem,length(eps1));
RCS_array_eps1_mnst_arr = zeros(colmn_elem,length(eps1));
N_arr = zeros(1,colmn_elem);
err = zeros(colmn_elem,length(eps1));
for j = 1:1:colmn_elem
    for i = 1:1:length(eps1)
        beam_array_eps1_mnst_arr(j,i) = sum(abs(Z_r(:,j))>eps1(i));
        RCS_array_eps1_mnst_arr(j,i) = 4*abs((sum(Z_r(abs(Z_r(:,j))>eps1(i),j))))^2 ;
    end
    err(j,:) = 20*log10(abs(1-RCS_array_eps1_mnst_arr(j,:)));
    err_temp = 20*log10(abs(1-RCS_array_eps1_mnst_arr(j,:)));
    if (err_temp(1)>-40)
        N = nan;
    else 
        N = beam_array_eps1_mnst_arr(j,find(err_temp>-40 == 1, 1));
        beam_array_eps1_mnst_arr(j,find(err_temp>-40 == 1, 1))
    end
    N_arr(j) = N;
end
%plot(beam_array_eps1_mnst_kmin8000,20*log10(abs(1-RCS_array_eps1_mnst_kmin8000)),'b','linewidth',2);
%hold on;
%plot(beam_array_eps1_bst_kmin8000,20*log10(abs(1-RCS_array_eps1_bst_kmin8000)),'r--','linewidth',2);
%grid on;
%xlabel('N');
%ylabel('Error(dB)');
%xlim([500,3000]);
%set(gca,'fontsize',20);
%xlim([200,5000]);
%%
%m = 6;
%plot(beam_array_eps1_mnst_arr(m,:),err_mnst(m,:),'linewidth',2);
%hold on 
%plot(beam_array_eps1_bst_arr(m,:),err_bst(m,:),'r--','linewidth',2);
%xlabel('{\it N}');
%ylabel('Error(dB)');
%xlim([0,5000]);
%ylim([-100 0]);
%set(gca,'fontsize',20);
%grid on;
