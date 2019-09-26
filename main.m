    %_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2011b(7.13)                                       %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%   Main paper: S. Mirjalili, A. Lewis                                    %
%               The Whale Optimization Algorithm,                         %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008   %
%                                                                         %
%_________________________________________________________________________%

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run WOA: [Best_score,Best_pos,WOA_cg_curve]=WOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________
clear all 
clc
clc;
% w=(0:0.01:pi)/pi;
w=linspace(0,pi,256);

%global fitness_wle ;

global j1 w beta;

% fitness_wle=20;
 
y=[0.05 -0.4 -1.1214 0.25]

 b_desired=y(1,1:2);
 a_desired=y(1,3:end);
 j1=freqz(b_desired,[1 a_desired],w);

dim=4;

% global fitness_wle count C Q P ws wp;
% count=0;
% fitness_wle=10;
% 
 SearchAgents_no=200; % Number of search agents
% 
% %Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% 
 Max_iteration=1000; % Maximum numbef of iterations
% dim=28;
% wp=0.4*pi;
% ws=0.6*pi;
lb=-1; ub=+1;
% 
% % Load details of the selected benchmark function
% %[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
%  N=dim;
% M=N/2;
% C=zeros(M,M);
% Q=zeros(M,M);
% P = zeros(M,1);
% F=zeros(M,M);
% for m=0:M-1
%     A=(N-1)/2-m;
%     for n=0:M-1
%         p=m-n;
%         q=m+n;
%         C(m+1,n+1)=ExpInt(p,q,N,ws,pi);
%         Q(m+1,n+1)=ExpInt(p,q,N,0,wp);
%     end;
%     P(m+1,1)=CosInt(A,0,wp);
% end;
for i=1:1
    tic;
[Best_score,Best_pos,WOA_cg_curve]=WOA(SearchAgents_no,Max_iteration,lb,ub,dim);
t(i,:)=toc
end

%  H(i,:)=[Best_score(i,1:M-1)/2 Best_score(i,M) fliplr(Best score(i,1:M-1)/2)];
  H=[Best_pos fliplr(Best_pos)]/2;
  H1= H/sum(H);
  [Hw w] =freqz(H1 ,1, 128);
% h_db=20*log10((abs(Hw)+eps)/max(abs(Hw)));
% delta_w=pi/length(Hw);
%Rp=-(min(h_db(1:wp/delta_w+1)))
%As=-(max(h_db(ws/delta_w+1:1:length(Hw))))
% % m=[1,1,1,0,0,0];
% % f=[0, Best_pos, 1];
% % [b,a]=yulewalk(6,f,m);
% % [h,w]=freqz(b,a,128)
% % plot(Best_pos,m,w/pi,abs(h))
fvtool(H1,1)
% [fmin minindex]=min(Best_score);
% result=[Best_score H]; 
% fvtool(H(minindex,:),1)
% save result_wle.mat Best_pos;
% save Fitness_wle.mat fitness_wle;



% figure('Position',[269   240   660   290])
% %Draw search space
% subplot(1,2,1);
% func_plot(Function_name);
% title('Parameter space')
% xlabel('x_1');
% ylabel('x_2');
% zlabel([Function_name,'( x_1 , x_2 )'])
% 
% %Draw objective space
% subplot(1,2,2);
semilogy(WOA_cg_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');
% 
% axis tight
% grid on
% box on
% legend('WOA')
% 
% display(['The best solution obtained by WOA is : ', num2str(Best_pos)]);
% display(['The best optimal value of the objective funciton found by WOA is : ', num2str(Best_score)]);

        



