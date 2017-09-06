% DugueSim.m
% quick simulation to check the analysis method of Dugue et al

nTrials = 1000;
N = 100;

%% Case where subjects perform at chance

% p1 and p2 are equal and performance is at .5
p1true = .5;
p2true = p1true;

P1 = nan(N,1);
P2 = nan(N,1);
Pdif = nan(N,1);

for n = 1:N
    
    resp = nan(nTrials,2);
    for i = 1:nTrials
        resp(i,1) = binornd(1,p1true);
        resp(i,2) = binornd(1,p2true);
    end
    resp = Shuffle(resp,1); % Shuffle so that these responses could have been on either side

    Pboth = sum(resp(:,1) == 1 & resp(:,2)== 1)/nTrials;
    Pnone = sum(resp(:,1) == 0 & resp(:,2)== 0)/nTrials;
    
    c = Pboth;
    b = 1 + Pboth - Pnone;
    
    det = b^2 - 4*c;
    
    P1(n) = (b + sign(det)*sqrt(abs(det)))/2;
    P2(n) = (b - sign(det)*sqrt(abs(det)))/2;
    
end

pc = mean(resp);
Pdif = P1-P2;

figure(1)
clf
hold on

subplot(4,3,1)
title('Divided attention')
plot(pc,'x')
xlim([0,3])
ylim([.25,1])
ylabel('Proportion correct')
xlabel('Left Right')

subplot(4,3,4)
hold on
hist(P1)
subplot(4,3,7)
hist(P2)

subplot(4,3,10)
hist(Pdif)
ylabel('Number of simulations')
xlabel('P1-P2')
xlim([-1,1])

%% Case where subjects perform equally in divided attention

% p1 and p2 are equal and performance is at .75
p1true = .75; 
p2true = p1true;

P1 = nan(N,1);
P2 = nan(N,1);
Pdif = nan(N,1);

for n = 1:N
    
    resp = nan(nTrials,2);
    for i = 1:nTrials
        resp(i,1) = binornd(1,p1true);
        resp(i,2) = binornd(1,p2true);
    end
    resp = Shuffle(resp,1); % Shuffle so that these responses could have been on either side
    
    Pboth = sum(resp(:,1) == 1 & resp(:,2)== 1)/nTrials;
    Pnone = sum(resp(:,1) == 0 & resp(:,2)== 0)/nTrials;
    
    b = 1 + Pboth - Pnone;
    c = Pboth;
    
    det = b^2 - 4*c;
    
    P1(n) = (b + sign(det)*sqrt(abs(det)))/2;
    P2(n) = (b - sign(det)*sqrt(abs(det)))/2;
    
end

pc = mean(resp);
Pdif = P1-P2;


subplot(4,3,2)
title('Divided attention')
plot(pc,'x')
xlim([0,3])
ylim([.25,1])
ylabel('Proportion correct')
xlabel('Left Right')

subplot(4,3,5)
hold on
hist(P1)
subplot(4,3,8)
hist(P2)

subplot(4,3,11)
hist(Pdif)
ylabel('Number of simulations')
xlabel('P1-P2')
xlim([-1,1])


%% Case where subjects perform unequally

% p1 and p2 are different and overall performance is .75
p1true = .5; 
p2true = .75;

P1 = nan(N,1);
P2 = nan(N,1);
Pdif = nan(N,1);

for n = 1:N
    
    resp = nan(nTrials,2);
    for i = 1:nTrials
        resp(i,1) = binornd(1,p1true);
        resp(i,2) = binornd(1,p2true);
    end
    resp = Shuffle(resp,1); % Shuffle so that these responses could have been on either side

    Pboth = sum(resp(:,1) == 1 & resp(:,2)== 1)/nTrials;
    Pnone = sum(resp(:,1) == 0 & resp(:,2)== 0)/nTrials;
    
    b = 1 + Pboth - Pnone;
    c = Pboth;
    
    det = b^2 - 4*c;
    
    P1(n) = (b + sign(det)*sqrt(abs(det)))/2;
    P2(n) = (b - sign(det)*sqrt(abs(det)))/2;
    
end

pc = mean(resp);
Pdif = P1-P2;


subplot(4,3,3)
title('Divided attention')
plot(pc,'x')
xlim([0,3])
ylim([.25,1])
ylabel('Proportion correct')
xlabel('Left Right')

subplot(4,3,6)
hold on
hist(P1)
subplot(4,3,9)
hist(P2)

subplot(4,3,12)
hist(Pdif)
ylabel('Number of simulations')
xlabel('P1-P2')
xlim([-1,1])