 %%% HW5 - BSS - Erfan Panahi 810198369
clc
clear
fprintf("HW5 - BSS - Erfan Panahi 810198369\n");

%% Definitions:

hw5 = load('hw5.mat');
TrainData_class1 = hw5.TrainData_class1;
TrainData_class1 = TrainData_class1 - mean(mean(TrainData_class1,3),2);
TrainData_class2 = hw5.TrainData_class2;
TrainData_class2 = TrainData_class2 - mean(mean(TrainData_class2,3),2);
TestData = hw5.TestData;
TestData = TestData - mean(mean(TestData,3),2);
TestLabel = hw5.TestLabel;

%% Part a.

[ChL,TimeL,N] = size(TrainData_class1);
RX1 = zeros(ChL); 
RX2 = zeros(ChL); 
for j = 1:N
    RX1 = RX1 + (TrainData_class1(:,:,j) * TrainData_class1(:,:,j).')./N;
    RX2 = RX2 + (TrainData_class2(:,:,j) * TrainData_class2(:,:,j).')./N;
end
[Wcsp , Dcsp] = eig(RX1,RX2);
Wcsp = Wcsp./(ones(length(Wcsp),1)*sqrt(sum(Wcsp.^2)));
Wcsp = Wcsp(:,end:-1:1);
Wcsp1 = Wcsp(:,1);
Wcsp30 = Wcsp(:,30);
X1 = TrainData_class1(:,:,49);
X2 = TrainData_class2(:,:,49);
WtX1_49_1 = Wcsp1.' * X1;
WtX2_49_1 = Wcsp1.' * X2;
WtX1_49_30 = Wcsp30.' * X1;
WtX2_49_30 = Wcsp30.' * X2;
t = 1:TimeL;
figure(1)
subplot(2,1,1);
plot(t,WtX1_49_1,'b',t,WtX2_49_1,'r');
title('$W_{csp}$ : var($W_1^TX$) for Test:49 (filter:$W_1$)','Interpreter','latex');
xlabel('time','Interpreter','latex');
ylabel('var($W_1^TX$)','Interpreter','latex');
legend('calss1','class2');
subplot(2,1,2);
plot(t,WtX1_49_30,'b',t,WtX2_49_30,'r');
title('$W_{csp}$ : var($W_{30}^TX$) for Test:49 (filter:$W_{30}$)','Interpreter','latex');
xlabel('time','Interpreter','latex');
ylabel('var($W_{30}^TX$)','Interpreter','latex');
legend('calss1','class2');
fprintf("\nPart a.\n  Test:49 , Filter:W1 ::: var(X1) = %f < var(X2) = %f\n",var(WtX1_49_1),var(WtX2_49_1));
fprintf("  Test:49 , Filter:W30 ::: var(X1) = %f > var(X2) = %f\n",var(WtX1_49_30),var(WtX2_49_30));

%% Part b.

figure(2)
subplot(2,1,1);
plot(1:ChL,abs(Wcsp1));
title('$W_{csp}$ : (filter:$W_1$)','Interpreter','latex');
xlabel('Channels','Interpreter','latex');
ylabel('$W_1$','Interpreter','latex');
subplot(2,1,2);
plot(1:ChL,abs(Wcsp30));
title('$W_{csp}$ : (filter:$W_{30}$)','Interpreter','latex');
xlabel('Channels','Interpreter','latex');
ylabel('$W_{30}$','Interpreter','latex');

%% Part c.

W_csp = Wcsp(:,[1:7,24:30]);
VarX1 = zeros(14,N);
VarX2 = zeros(14,N);
for i = 1:N
    x1 = TrainData_class1(:,:,i);
    x2 = TrainData_class2(:,:,i);
    WtX1 = W_csp.' * x1;
    WtX2 = W_csp.' * x2;
    VarX1(:,i) = (var(WtX1.')).';
    VarX2(:,i) = (var(WtX2.')).';
end
mu1 = mean(VarX1,2);
mu2 = mean(VarX2,2);
Sigma1 = zeros(14);
Sigma2 = zeros(14);
for i = 1:N
    Sigma1 = Sigma1 + ((VarX1(:,i)-mu1)*(VarX1(:,i)-mu1).')./N;
    Sigma2 = Sigma2 + ((VarX2(:,i)-mu2)*(VarX2(:,i)-mu2).')./N;
end
A = (mu1-mu2)*(mu1-mu2).';
B = Sigma1 + Sigma2;
[Wlda,Dlda] = eig(A,B);
Wlda = Wlda(:,end:-1:1);
W_lda = Wlda(:,1);
mu_1 = W_lda.' * mu1;
mu_2 = W_lda.' * mu2;
c = (mu_1 + mu_2)/2;
fprintf("\nPart c.\n WldaT = %s\n c = %f\n",num2str(W_lda.'),c);


%% Part d.

fprintf("\nPart d.\n mu2 = %f < c = %f < mu1 = %f\n\n",mu_2,c,mu_1);
[C,T,NumTest] = size(TestData);
Var_Test = zeros(14,NumTest);
for i = 1:NumTest
   WtX_Test = W_csp.' * TestData(:,:,i);
   Var_Test(:,i) = var(WtX_Test.').';
end
mu = W_lda.' * Var_Test;
Label = 2*(mu < c) + (mu > c);

%% Part e.

figure(3)
plot(1:NumTest,Label,'rx',1:NumTest,TestLabel,'bo');
ylim([0.5,2.5]);
title('Estimated and Actual Labels','Interpreter','latex');
xlabel('Number of Test','Interpreter','latex');
ylabel('Label','Interpreter','latex');
legend('Estimated','Actual');
