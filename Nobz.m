%% Number of frequency of the beams as a function of frequency and the monostatic and the bistatic RCS z_0 = -2a
%% Monostatic RCS z_min = -a 
semilogy(log2(k_sub/1000), N_arr_mnst_zmina_xi,'-+m','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_mnst_zmina_xix,'-xk','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_mnst_zmina_xxi,'-ob','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_mnst_zmina_x,'-*r','linewidth',2);
grid on
p_1 = semilogy(log2(k_sub(5:8)/1000), N_arr_mnst_zmina_xxi(5:8),'ob','linewidth',2);
p_1.MarkerSize = 10;
p_1 = semilogy(log2(k_sub(5:8)/1000), N_arr_mnst_zmina_xi(5:8),'+m','linewidth',2);
p_2 = semilogy(log2(k_sub(5:8)/1000), N_arr_mnst_zmina_xi(5:8),'+m','linewidth',2);
p_2.MarkerSize = 10;
p_3 = semilogy(log2(k_sub(9:12)/1000), N_arr_mnst_zmina_xix(9:12),'xk','linewidth',2);
p_3.MarkerSize = 10;
xticks([0,1,2,3,4])
xlabel('{\it j}');
ylabel('{\it N}');
xlim([0,4])
ylim([100 10^6]);
l_1 = legend('$\xi$','$\xi$-$\textit{x}$','$\textit{x}$-$\xi$','$\textit{x}$');
set(l_1,'Interpreter','latex');
set(gca,'fontsize',20);
%% Bistatic RCS z_min = -a 
semilogy(log2(k_sub/1000), N_arr_bst_zmina_xi,'-+m','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_bst_zmina_xix,'-xk','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_bst_zmina_xxi,'-ob','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_bst_zmina_x,'-*r','linewidth',2);
grid on
p_1 = semilogy(log2(k_sub(5:8)/1000), N_arr_bst_zmina_xxi(5:8),'ob','linewidth',2);
p_1.MarkerSize = 10;
p_1 = semilogy(log2(k_sub(5:8)/1000), N_arr_bst_zmina_xi(5:8),'+m','linewidth',2);
p_2 = semilogy(log2(k_sub(5:8)/1000), N_arr_bst_zmina_xi(5:8),'+m','linewidth',2);
p_2.MarkerSize = 10;
p_3 = semilogy(log2(k_sub(9:12)/1000), N_arr_bst_zmina_xix(9:12),'xk','linewidth',2);
p_3.MarkerSize = 10;
xticks([0,1,2,3,4])
xlabel('{\it j}');
ylabel('{\it N}');
xlim([0,4])
ylim([100 10^6]);
l_1 = legend('$\xi$','$\xi$-$\textit{x}$','$\textit{x}$-$\xi$','$\textit{x}$');
set(l_1,'Interpreter','latex');
%legend('\xi','\xi-x','{\it x}-\xi','x')
set(gca,'fontsize',20);
%% Bistatic RCS z_0 = 0
semilogy(log2(k_sub/1000), N_arr_bst_z0_xi,'-+m','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_bst_z0_xix,'-xk','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_bst_z0_xxi,'-ob','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_bst_z0_x,'-*r','linewidth',2);
grid on
p_1 = semilogy(log2(k_sub(5:8)/1000), N_arr_bst_z0_xxi(5:8),'ob','linewidth',2);
p_1.MarkerSize = 10;
p_1 = semilogy(log2(k_sub(5:8)/1000), N_arr_bst_z0_xi(5:8),'+m','linewidth',2);
p_2 = semilogy(log2(k_sub(5:8)/1000), N_arr_bst_z0_xi(5:8),'+m','linewidth',2);
p_2.MarkerSize = 10;
p_3 = semilogy(log2(k_sub(9:12)/1000), N_arr_bst_z0_xix(9:12),'xk','linewidth',2);
p_3.MarkerSize = 10;
xticks([0,1,2,3,4])
xlabel('{\it j}');
ylabel('{\it N}');
xlim([0,4])
ylim([100 10^6]);
l_1 = legend('$\xi$','$\xi$-$\textit{x}$','$\textit{x}$-$\xi$','$\textit{x}$');
set(l_1,'Interpreter','latex');
%legend('\xi','\xi-x','{\it x}-\xi','x')
set(gca,'fontsize',20);
xticks([0,1,2,3,4])
xlabel('{\it j}');
ylabel('{\it N}');
xlim([0,4])
ylim([100 10^6]);
l_1 = legend('$\xi$','$\xi$-$\textit{x}$','$\textit{x}$-$\xi$','$\textit{x}$');
set(l_1,'Interpreter','latex');
%legend('\xi','\xi-x','{\it x}-\xi','x')
set(gca,'fontsize',20);
%% Monostatic RCS z_0 = 0 
semilogy(log2(k_sub/1000), N_arr_mnst_z0_xi,'-+m','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_mnst_z0_xix,'-xk','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_mnst_z0_xxi,'-ob','linewidth',2);
hold on
semilogy(log2(k_sub/1000), N_arr_mnst_z0_x,'-*r','linewidth',2);
grid on
p_1 = semilogy(log2(k_sub(5:8)/1000), N_arr_mnst_z0_xxi(5:8),'ob','linewidth',2);
p_1.MarkerSize = 10;
p_1 = semilogy(log2(k_sub(5:8)/1000), N_arr_mnst_z0_xi(5:8),'xm','linewidth',2);
p_2 = semilogy(log2(k_sub(5:8)/1000), N_arr_mnst_z0_xi(5:8),'xm','linewidth',2);
p_2.MarkerSize = 10;
p_3 = semilogy(log2(k_sub(9:12)/1000), N_arr_bst_z0_xix(9:12),'xk','linewidth',2);
p_3.MarkerSize = 10;
xticks([0,1,2,3,4])
xlabel('{\it j}');
ylabel('{\it N}');
xlim([0,4])
ylim([100 10^6]);
l_1 = legend('$\xi$','$\xi$-$\textit{x}$','$\textit{x}$-$\xi$','$\textit{x}$');
set(l_1,'Interpreter','latex');
%legend('\xi','\xi-x','{\it x}-\xi','x')
set(gca,'fontsize',20);




