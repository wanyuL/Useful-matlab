%clear; clc;


%% SIMULATE THE PROCESS USING THE GILLESPIE ALGORITHM.
%Modeling single-gene expression. mRNA molecules are transcribed at rate
%kR from the template DNA strand. Proteins are translated at a rate kP from
%each mRNA molecule. Proteins and mRNA degrade at rates gP and gR,
%respectively. 

%define our parameters
gR=0.1;     %1/sec
gP=0.002;   %1/sec
kR=0.01;    %1/sec
kP=1;       %1/sec  i.e. b=kP/gR

Khalf = 16;
hill = 5;
kLeak = 0.001;


%define the reaction propensities
Rsyn = @(r,p) kR*p^hill/(p^hill+Khalf^hill)+kLeak;
Rdeg = @(r,p) gR*r;
Psyn = @(r,p) kP*r;
Pdeg = @(r,p) gP*p;

%define the time variable and the end time
time = 0;
endTime = 3600 * 100; % in seconds

%initial condition
mRNA = 0;
protein = 0;

% define when we should output, and a matrix to store the results
dtForOutput = 1; % output every second
nextOutputTime = time+dtForOutput; % set the time to output next
results = zeros(endTime/dtForOutput,3); %array to store output
results(1,:) = [time mRNA protein];
lastResults = results(1,:);
resultIndex = 2;

%define the list of reaction propensities based on the initial
%conditions and parameters.
propensities = [ ...
    Rsyn(mRNA,protein)... %reaction 1: transcription of a gene into a transcript
    Rdeg(mRNA,protein)... %reaction 2: degradation of a transcript
    Psyn(mRNA,protein)... %reaction 3: translation of a transcript into a protein
    Pdeg(mRNA,protein)... %reaction 4: degradation of a protein
    ];

% generate our random number ahead of time for efficiency, save it
% to a matrix called randNums.  We also need a counter to determine
% which random number to use: Gillespie's algorithm uses 2 random numbers
% at each step.
numberOfRandomNums = 5000000;
randNums = rand(numberOfRandomNums,2); 
rIndex = 1;

% begin the main simulation loop of the algorithm
while time < endTime
    
    % step 0, if we used up all random numbers, generate a new batch
    if ( rIndex > numberOfRandomNums )
        fprintf('time = %4.2f Used all random numbers. Generating new batch\n',time);
        randNums = rand(numberOfRandomNums,2); 
        rIndex = 1;
    end
    
    % step 1, generate the time to the next event (requires us
    % to calculate the total propensity of the system! The total 
    % propensity is used to generate a random number exponentialy 
    % distributed). 
    rtot = sum(propensities);
    tau = -(1/rtot)*log(randNums(rIndex,1));
    
    % step 2, select the next reaction to fire- first we have to
    % draw a random number between 0 and the total propensity sum
    %value = rand(1)*rtot; %
    value = randNums(rIndex,2)*rtot;
    
    % step2 / step 3, actually fire the reaction by updating the
    % counts based on 'value' that we generated.  
    
    % based on the value and the running sum, we can fire the selected rxn.
    if ( value < propensities(1) )
        % reaction 1 fires: transcription
        mRNA = mRNA + 1;
       % fprintf('r1\n');
    elseif ( value < propensities(1)+propensities(2) )
        % reaction 2: degradation of one transcript
        mRNA = mRNA - 1;
       % fprintf('r2\n');
    elseif ( value < propensities(1)+propensities(2)+propensities(3) )
        % reaction 3: translation
        protein = protein + 1;
       % fprintf('r2\n');
    elseif ( value < propensities(1)+propensities(2)+propensities(3)+propensities(4) )
        % reaction 4: degradation of a protein
        protein = protein - 1;
       % fprintf('r2\n');
    else
        % we should never get here, but we may if there are numerical
        % errors or if you calculated your propensities incorrectly
        fprintf('numerical errors.  quitting'); break;  
    end


    % step 4, update the propensities and the simulation time
    propensities = [ 
        Rsyn(mRNA,protein) ...
        Rdeg(mRNA,protein) ...
        Psyn(mRNA,protein) ...
        Pdeg(mRNA,protein) ];
    time = time + tau;
    rIndex = rIndex+1;
    
    % step 5, generate results from the simulation by saving it to 
    % our results output matrix.  Only generate output when the simulation
    % time exceeds the nextOutputTime.  
    if time>nextOutputTime
        
        % if we jumped over the nextOutputTime we need to update our array
        while time>nextOutputTime
           results(resultIndex,:) = [nextOutputTime lastResults(2) lastResults(3)];
           nextOutputTime = nextOutputTime + dtForOutput;
           resultIndex = resultIndex+1; 
           if nextOutputTime>endTime
               break;
           end;
        end;
        
    end
        
    % store results of last jump in case next we jump beyond storage point 
    % and need to store intermediate time points
    lastResults = [time mRNA protein];
    
end;

%% PLOTS
Rmax = 8;   %max(results(:,2))+1;
Pmax = 162; %max(results(:,3))+1;
tmax = 3600*50;

tt=results(:,1); %time axis

figure;
linestyle = 'k-';
subplot(2,2,1); 
plot(tt/60,results(:,2),linestyle);
xlabel('time(min)'); ylabel('mRNA (#)');
axis([0,tmax/60,0,Rmax]);
title('Time dependent dynamics');
hold on; 
plot(tt/60,(kR/gR)*(1-exp(-gR*tt)),'r');
hold off;

subplot(2,2,3); 
plot(tt/60,results(:,3),linestyle);
xlabel('time(min)'); ylabel('protein (#)');
axis([0,tmax/60,0,Pmax]);
title('Time dependent dynamics');
hold on; 
plot(tt/60,(kR*kP/gR/gP).*(1-(gP.*exp(-gR.*tt)-gR.*exp(-gP.*tt))./(gP-gR)),'r');
hold off;

%If we wait long enough (How much???) the system will forget about the
%initial condition and the system will tend towards the steady state
%distribution. 

ndrop=3600*10; %number of time points to drop to ensure system reached steady state 
Rcv = sqrt(var(results(ndrop:end,2)))/mean(results(ndrop:end,2));
Pcv = sqrt(var(results(ndrop:end,3)))/mean(results(ndrop:end,3));

subplot(2,2,2); 
[pR,Rbins] = hist(results(ndrop:end,2),0:max(results(:,2)));
pR = pR/sum(pR);
barh(Rbins,pR);
ylabel('mRNA (#)'); xlabel('Probability');
ylim([0,Rmax]);
title(strcat('Steady state distribution. CV=',num2str(Rcv,2)));

subplot(2,2,4); 
[pP,Pbins] = hist(results(ndrop:end,3),0:max(results(:,3)));
pP = pP/sum(pP);
barh(Pbins,pP);
ylabel('protein (#)'); xlabel('Probability');
ylim([0,Pmax]);
title(strcat('Steady state distribution. CV=',num2str(Pcv,2)));
