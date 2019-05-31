%% Figure 3
%% Savings experiment (with two guassian input profiles)
%% (using Brian simulations)

clc;
clear all;
close all;

%global dt;
LTMSD_fonts;
fontsize = 11;

save_on = 0;
save_path = '';
save_filename = ['Fig3savings'];

%% Load results from Brian output
[brianParams] = importdata([save_path 'fromBrian/outParams.br']);

N = 1;

resolution_export = brianParams(2);
nruns = brianParams(end-1);

aux = importdata([save_path 'fromBrian/outA_run0.br']);

brianA = zeros(nruns, size(aux,1), size(aux,2));
brianU = zeros(nruns, size(aux,1), size(aux,2));
for i=1:nruns
    [brianA(i,:,:)] = importdata([save_path 'fromBrian/outA_run' num2str(i-1) '.br']);
    [brianU(i,:,:)] = importdata([save_path 'fromBrian/outU_run' num2str(i-1) '.br']);
end

brianAmean = reshape(mean(brianA,1), size(brianA,2),size(brianA,3));
brianAsem = reshape(std(brianA,0,1)./sqrt(nruns), size(brianA,2),size(brianA,3));
brianUmean = reshape(mean(brianU,1), size(brianA,2),size(brianA,3));
brianUsem = reshape(std(brianU,0,1)./sqrt(nruns), size(brianA,2),size(brianA,3));

%brianA = brianAmean;
%brianU = brianUmean;

%% Figure
figure('Position',[10 250 1200 450]);
set(gcf, 'Units', 'centimeters');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
a4 = get(gcf, 'PaperSize');
set(gcf, 'Position', [1 10 a4(1) a4(2)-10]);
set(gcf, 'PaperSize', [a4(1) a4(2)-10]);


    nrows = 4;
    ncols = 4;
    % Matlab figure representation
    fig_ind = reshape(1:nrows*ncols, ncols, nrows)';
    A = 1;
    B = A+1;
    C = B+1;
    D = C+1;
    E = D+1;
    
    fig_subind = [A A A A;
                  B B B B;
                  C C C C;
                  D D E E];



%% A - Toy simulation

    %Options
    %n_shifts = 4;
    shift = 10000; %ms
    shift_short = 1000; %ms
    %shift = 200; %ms
    %stime = shift*dt * n_shifts; %Simulation time (s) - Monocular deprivation (eye 1: low, eye2: high)
    
   % stime2 = 3; %Simulation time (s) - Binocolar vision (eye 1: high, eye 2: low)
   % stime3 = 3; %Simulation time (s) - Monocular deprivation (eye 1: low, eye2: high) 
    n_conditions = 1;
    resolution = 1;
    verbose = 1;
    
    %rates = [5 50]; %Default rate for the poisson input
    %rates2 = [50 5];
    
    Ainit = 1; %2e-3; %1 for network simulations
    Amin = 0.1;
    Amax = brianParams(end);
    decay = 0.1;
    
    Uinit = 0.65;
    Umin = 0.1;
    Umax = 1;
    lr = 0.25;
    
    n_syn = brianParams(1); %N synapses
    n_synI = 0; %Ratio of interneurons
    
    %AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp tetaFB tetaFF
    p_on = [1 1 1 1 1 1 0 0];

    %%Short-term Plasticity
    %(tau_r, tau_p, tau_in, Pb, tau_mem, facil, Ase)
    tau_in = 2e-3;
    tau_mem = 20e-3;
    
    if(exist('Ase', 'var')==0)
        Ase = 5e-9; %Synaptic efficacy
        %Ase = 9e-10; %Synaptic efficacy
        %Ase = 0; %Synaptic efficacy
    end    
    %Ase = 25.0e-7; %Synaptic efficacy
    AseD = Ase; %5e-9;
    AseF = Ase; %5e-9;
    
         
    STD_params = [500e-03, 5e-3, Umin, 0.06]; %Depression
    
    %N gaussian profiles
    s = 5;
    ngauss = 2;
    
    if(ngauss>1)
        pos = [brianParams(5) brianParams(6)];
    else
        pos = n_syn/2;
    end
    rmin = 5;
    rmax = 50;

    x = 1:n_syn;
    rates = zeros(ngauss,length(x));
    Ui = ones(1,length(x)).*Umin;
    Ai = ones(1,length(x));

    stime = brianParams(3)*1e3;         %Simulation time (s)
    stime2 = brianParams(4)*1e3;         %Simulation time (s)
    shiftf = [stime stime2 stime2 stime2];
    
    var_bnoise = 0.35;
    
    %9. Plot 
    c = hot(n_conditions*n_syn*(1-n_synI)+2);
    c = [[0 0 0]; c];
    c2 = 1;
    leg = {};
    VARa = zeros(n_syn, size(brianAmean, 2));
    SNRa = zeros(n_syn, size(brianAmean, 2));
    Wa = zeros(n_syn, size(brianAmean, 2));
    Aa = zeros(n_syn, size(brianAmean, 2));
    Ua = zeros(n_syn, size(brianAmean, 2));
    
%     postRates = zeros(n_syn, stime/dt-1);
%     g = 0.1;
%     lE = 50;
%     aI = 0;
%     lI = 50;
%     tau_r = 2e-3;
%     tau_m = 20e-3;
%     Vreset = 0;%0;
%     Vth = 20e-3;%10e-3;
%     u = Vreset:0.01e-3:Vth;
    
    for j=1:n_syn
        if(j>(n_synI*n_syn))
            VARa(j,:) = N.*s.*(brianAmean(j,:).^2).*(1-brianUmean(j,:));
            SNRa(j,:) = 2*((N.*brianUmean(j,:).*brianAmean(j,:).^2)./(N.*brianUmean(j,:).*(brianAmean(j,:).^2).*(1-brianUmean(j,:))+var_bnoise));
            %Wa(j,:) = prules{i,j}.wA(1:end-1);
            Aa(j,:) = brianAmean(j,:);
            Ua(j,:) = brianUmean(j,:);

            %Calc post firing rates using analytical results
%                 for t=1:stime/dt-1
%                     aE = -g*prules{i,j}.PA(t)*prules{i,j}.UA(t);
% 
%                     muQ = tau_m * (aE*lE + aI*lI);
%                     sigmasQ = (tau_m/2) * (aE^2*lE + aI^2*lI);
% 
%                     intu = trapz(u, exp((u-muQ).^2./(2*sigmasQ)) .* (1 + erf((u-muQ)./(sqrt(sigmasQ)*sqrt(2)))));
% 
%                     postRates(j, t) = 1/(tau_r + (tau_m/sqrt(sigmasQ)) * sqrt(pi/2) * intu);
%                 end

            c2=c2+1;
        end
    end

    
    %Panel A
    fig_a = subplot(nrows, ncols, fig_ind(fig_subind==A), 'replace');
    p2 = get(fig_a, 'Position');
    p2(2) = p2(2)+0.04;
    set(fig_a, 'Position', p2);
    
    imagesc(1:size(Ua,2), 1:size(Ua,1), Ua);
    hold on;
    
    p = 75;
    %scatter(fig_a, shift_short, p, 120, [1 0.5 0.5], 'filled');
    %scatter(fig_a, (shift_short+shift+shift), p, 120, [198/255 17/255 0], 'filled');
    
    set(fig_a, 'XTickLabel', []);
    cbh = colorbar;
    set(get(cbh,'title'),'string','P','fontsize',fontsize);
    caxis([Umin Umax]);
    xlim(fig_a, [0 size(Ua,2)-100]);
    yLabA = ylabel(fig_a, 'input number', 'FontSize', fontsize);
    pyLabA = get(yLabA, 'Position');
    set(yLabA, 'Position', [pyLabA(1)+125 pyLabA(2) pyLabA(3)]);
    set(fig_a, 'TickDir', 'Out');
    
    %Panel B
    fig_b = subplot(nrows, ncols, fig_ind(fig_subind==B), 'replace');
    p2 = get(fig_b, 'Position');
    p2(2) = p2(2)+0.06;
    set(fig_b, 'Position', p2);
    
    imagesc(1:size(Aa,2),1:size(Aa,1), Aa);
    %xlabel(fig_b, 'Time (ms)', 'FontSize', fontsize);
    set(fig_b, 'XTickLabel', []);
    ylabel(fig_b, 'input number', 'FontSize', fontsize);
    cbh = colorbar;
    set(get(cbh,'title'), 'string', 'q', 'fontsize', fontsize);
    caxis([Amin Amax]);
    
    colormap(flipdim(gray,1));
    xlim(fig_b, [0 size(Aa,2)-100]);
    set(fig_b, 'TickDir', 'Out');
    %keyboard

shiftpre = shift;
    win = 100;
    
    %% C-D - Learning curve (Postsynaptic freq)
    
    fig_c = subplot(nrows, ncols, fig_ind(fig_subind==C), 'replace');
        p2 = get(fig_c, 'Position');
        p2(2) = p2(2)+0.075;
        p2(3) = p2(3)-0.087;
        set(fig_c, 'Position', p2);
        
        aux = importdata([save_path 'fromBrian/rateInput1_run0.br']);

        brianRateInput1 = zeros(nruns, size(aux,1), size(aux,2));
        brianRateInput2 = zeros(nruns, size(aux,1), size(aux,2));
        for i=1:nruns
            [brianRateInput1(i,:,:)] = importdata([save_path 'fromBrian/rateInput1_run' num2str(i-1) '.br']);
            [brianRateInput2(i,:,:)] = importdata([save_path 'fromBrian/rateInput2_run' num2str(i-1) '.br']);
        end
        
        brianInput1mean = reshape(mean(brianRateInput1,1), size(brianRateInput1,2),size(brianRateInput1,3));
        brianInput1sem = reshape(std(brianRateInput1,0,1)./sqrt(nruns), size(brianRateInput1,2),size(brianRateInput1,3));
        brianInput2mean = reshape(mean(brianRateInput2,1), size(brianRateInput2,2),size(brianRateInput2,3));
        brianInput2sem = reshape(std(brianRateInput2,0,1)./sqrt(nruns), size(brianRateInput2,2),size(brianRateInput2,3));

        brianRateInput1 = brianInput1mean;
        brianRateInput2 = brianInput2mean;

        xx = (1:length(brianRateInput1)).*win/1e3;
        
         %figure
         %x=0:0.01:2*pi;                  %#initialize x array
         %y1=sin(x);                      %#create first curve
         %y2=sin(x)+.5;                   %#create second curve
         %X=[x,fliplr(x)];                %#create continuous x value array for plotting
         %Y=[y1,fliplr(y2)];              %#create y values for out and then back
         %fill(X,Y,'b');                  %#plot filled area
        
        y1=max(brianInput1mean-brianInput1sem,0)';
        y2=max(brianInput1mean+brianInput1sem,0)';
        X=[xx, fliplr(xx)];
        Y=[y1, fliplr(y2)];
        h=fill(X,Y, red2, 'EdgeColor', [1 1 1]);
        set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
        hold on;
        
        %h=area(fig_c, xx', max([brianInput1mean-brianInput1sem 2.*brianInput1sem],0));
        %set(h(2),'FaceColor', red2,'EdgeColor', [1 1 1]);
        %set(h(1),'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
        %set(get(get(h(1),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
        %set(get(get(h(2),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
        %hold on
%         h=area(fig_c, xx', max([brianInput2mean-brianInput2sem 2.*brianInput2sem], 0));
%         set(h(2),'FaceColor', blue2, 'EdgeColor', [1 1 1]);
%         set(h(1),'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
%         set(get(get(h(1),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
%         set(get(get(h(2),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
        
        y1=max(brianInput2mean-brianInput2sem,0)';
        y2=max(brianInput2mean+brianInput2sem,0)';
        Y=[y1, fliplr(y2)];
        h=fill(X, Y, blue2, 'EdgeColor', [1 1 1]);
        set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend

        plot(fig_c, xx, brianInput1mean, '-', 'Color', red1, 'LineWidth', lineWidth);
        plot(fig_c, xx, brianInput2mean', '-', 'Color', blue1, 'LineWidth', lineWidth);
        
        set(fig_c,'layer','top')
        xlabel(fig_c, 'time (s)', 'FontSize', fontsize);
        ylabel(fig_c, 'post rate (Hz)', 'FontSize', fontsize);
        
        legC = legend(fig_c, 'Stimulus 1', 'Stimulus 2', 'Location', 'SouthEast');
        legend(fig_c, 'boxoff');
        posC = get(legC, 'Position');
        set(legC, 'Position', [posC(1)+0.075 posC(2)+0.02 posC(3) posC(4)]);
        
        box off;
        ylim(fig_c, [0 50.5]);
        xlim(fig_c, [0 xx(end)-2.5]);
        set(fig_c, 'TickDir', 'Out');
        
    if(0)
        
    fig_d = subplot(nrows, ncols, fig_ind(fig_subind==D), 'replace');
        pd = get(fig_d, 'Position');
        pd(3) = pd(3)-0.03;
        set(fig_d, 'Position', pd);
        
        
        xx = (0:win:stime2)/1e3;
        
        y1=max(brianInput2mean(stime/win:(stime+stime2)/win)-brianInput2sem(stime/win:(stime+stime2)/win),0)';
        y2=max(brianInput2mean(stime/win:(stime+stime2)/win)+brianInput2sem(stime/win:(stime+stime2)/win),0)';
        X=[xx, fliplr(xx)];
        Y=[y1, fliplr(y2)];
        h=fill(X,Y, blue2, 'EdgeColor', [1 1 1]);
        set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
        
        hold on
        
        y1=max(brianInput2mean((stime+stime2*2)/win:(stime+stime2*3)/win)-brianInput2sem((stime+stime2*2)/win:(stime+stime2*3)/win),0)';
        y2=max(brianInput2mean((stime+stime2*2)/win:(stime+stime2*3)/win)+brianInput2sem((stime+stime2*2)/win:(stime+stime2*3)/win),0)';
        
        Y=[y1, fliplr(y2)];
        h=fill(X,Y, black2, 'EdgeColor', [1 1 1]);
        set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
        
        
        
%         h=area(fig_d, xx', max([brianInput2mean((stime+stime2*2)/win:(stime+stime2*3)/win)-brianInput2sem((stime+stime2*2)/win:(stime+stime2*3)/win) 2.*brianInput2sem((stime+stime2*2)/win:(stime+stime2*3)/win)], 0));
%         set(h(2),'FaceColor', black2, 'EdgeColor', [1 1 1]);
%         set(h(1),'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
%         set(get(get(h(1),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
%         set(get(get(h(2),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
        
        %hold on
        
%         h=area(fig_d, xx', max([brianInput2mean(stime/win:(stime+stime2)/win)-brianInput2sem(stime/win:(stime+stime2)/win) 2.*brianInput2sem(stime/win:(stime+stime2)/win)],0));
%         set(h(2),'FaceColor', red2,'EdgeColor', [1 1 1]);
%         set(h(1),'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
%         set(get(get(h(1),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
%         set(get(get(h(2),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
        
    
        %Learning
        plot(fig_d, xx, brianInput2mean(stime/win:(stime+stime2)/win), 'Color', blue1, 'LineWidth', lineWidth);
        %Relearning
        plot(fig_d, xx, brianInput2mean((stime+stime2*2)/win:(stime+stime2*3)/win), 'Color', black1, 'LineWidth', lineWidth);

        set(fig_d,'layer','top')
        xlabel(fig_d, 'time (s)', 'FontSize', fontsize);
        ylabel(fig_d, 'post. rate (Hz)', 'FontSize', fontsize);
        legD = legend(fig_d, 'Learning', 'Relearning', 'Location', 'NorthEast');
        posD = get(legD, 'Position');
        posD(2) = posD(2)+0.02;
        set(legD, 'Position', posD);
        legend(fig_d, 'boxoff');
        
        box off;
        ylim(fig_d, [0 85]);
        xlim(fig_d, [0 xx(end)-2.5]);
        set(fig_d, 'TickDir', 'Out');
    end
    
if(exist('shiftpre', 'var'))
    shift = shiftpre;    
else
    shift = stime;
end



%% D - Learning curve (P(discriminability))
fig_d = subplot(nrows, ncols, fig_ind(fig_subind==D), 'replace');
        pd = get(fig_d, 'Position');
        pd(3) = pd(3)-0.01;
        pd(1) = pd(1)+0.02;
        pd(2) = pd(2)+0.03;
        set(fig_d, 'Position', pd);

%Calculate P(discriminability) in p

var_bnoise = 0.35;
T = -3:0.1:3; %Dissimilarity threshold

p_fp = 1/2*erfc(T./sqrt(2*var_bnoise));
p_fp(isnan(p_fp)) = 0;


stime2r = stime2/resolution_export;
stimer = stime/resolution_export;


%1st switch
p_tp = zeros(nruns, stime2r, size(T,2));
p_tp_fixedU = zeros(nruns, stime2r, size(T,2));
for i=1:nruns
    for t=1:stime2r
        p_tp(i,t,:) = 1/2*erfc((T-(N*brianU(i, p, stimer+t)*brianA(i, p, stimer+t)))./sqrt(2*N*brianU(i,p, stimer+t)*brianA(i, p, stimer+t)^2*(1-brianU(i, p, stimer+t))));    
        p_tp_fixedU(i,t,:) = 1/2*erfc((T-(N*brianU(i, p, 1)*brianA(i, p, stimer+t)))./sqrt(2*N*brianU(i, p, 1)*brianA(i, p, stimer+t)^2*(1-brianU(i,p, 1))));    
    
        p_tp(i, t,isnan(p_tp(i,t,:))) = 0;
        p_tp_fixedU(i, t,isnan(p_tp_fixedU(i,t,:))) = 0;
    end
end    


res_lstp_1st_switch = zeros(nruns, stime2r,1);
res_lstp_1st_switch_fixedU = zeros(nruns, stime2r,1);
for i=1:nruns
    for t=1:stime2r
        res_lstp_1st_switch(i,t) = trapz(flipdim(p_fp,2), flipdim(p_tp(i,t,:),2));
        res_lstp_1st_switch_fixedU(i,t) = trapz(flipdim(p_fp,2), flipdim(p_tp_fixedU(i,t,:),2));
    end
end

res_lstp_1st_switch_mean = reshape(mean(res_lstp_1st_switch,1), size(res_lstp_1st_switch,2),size(res_lstp_1st_switch,3));
res_lstp_1st_switch_sem = reshape(std(res_lstp_1st_switch,0,1)./sqrt(nruns), size(res_lstp_1st_switch,2),size(res_lstp_1st_switch,3));

xx = (1:resolution_export:stime2)./1e3;

%2nd switch
p_tp = zeros(nruns, stime2r, size(T,2));
p_tp_fixedU = zeros(nruns, stime2r, size(T,2));
for i=1:nruns
    for t=1:stime2r
        p_tp(i, t,:) = 1/2*erfc((T-(N*brianU(i,p, (stimer+stime2r+stime2r)+t)*brianA(i, p, (stimer+stime2r+stime2r)+t)))./sqrt(2*N*brianU(i, p, (stimer+stime2r+stime2r)+t)*brianA(i, p, (stimer+stime2r+stime2r)+t)^2*(1-brianU(i, p, (stimer+stime2r+stime2r)+t))));    
        p_tp_fixedU(i, t,:) = 1/2*erfc((T-(N*brianU(i, p, 1)*brianA(i, p, (stimer+stime2r+stime2r)+t)))./sqrt(2*N*brianU(i, p, 1)*brianA(i, p, (stimer+stime2r+stime2r)+t)^2*(1-brianU(i, p, 1))));    
    
        p_tp(i, t, isnan(p_tp(i, t,:))) = 0;
        p_tp_fixedU(i, t, isnan(p_tp_fixedU(i, t,:))) = 0;
    end
end    

res_lstp_2nd_switch = zeros(nruns, stime2r,1);
res_lstp_2nd_switch_fixedU = zeros(nruns, stime2r,1);
for i=1:nruns
    for t=1:stime2r
        res_lstp_2nd_switch(i,t) = trapz(flipdim(p_fp,2), flipdim(p_tp(i,t,:),2));
        res_lstp_2nd_switch_fixedU(i,t) = trapz(flipdim(p_fp,2), flipdim(p_tp_fixedU(i,t,:),2));
    end
end    

res_lstp_2nd_switch_mean = reshape(mean(res_lstp_2nd_switch,1), size(res_lstp_2nd_switch,2),size(res_lstp_2nd_switch,3));
res_lstp_2nd_switch_sem = reshape(std(res_lstp_2nd_switch,0,1)./sqrt(nruns), size(res_lstp_2nd_switch,2),size(res_lstp_2nd_switch,3));

y1=max(res_lstp_2nd_switch_mean-res_lstp_2nd_switch_sem,0)';
y2=max(res_lstp_2nd_switch_mean+res_lstp_2nd_switch_sem,0)';
X=[xx, fliplr(xx)];
Y=[y1, fliplr(y2)];

%Random choice
h=plot(fig_d, [0 xx(end)-9], [0.5 0.5], ':k', 'Color', black1, 'LineWidth', 1);
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
hold on;

h=fill(X, Y, black2, 'EdgeColor', [1 1 1]);
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
hold on;

y1=max(res_lstp_1st_switch_mean-res_lstp_1st_switch_sem,0)';
y2=max(res_lstp_1st_switch_mean+res_lstp_1st_switch_sem,0)';
Y=[y1, fliplr(y2)];
h=fill(X, Y, blue2, 'EdgeColor', [1 1 1]);
set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
      
% h=area(fig_e, xx', max([res_lstp_1st_switch_mean-res_lstp_1st_switch_sem 2.*res_lstp_1st_switch_sem],0));
% set(h(2),'FaceColor', red2,'EdgeColor', [1 1 1]);
% set(h(1),'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
% set(get(get(h(1),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
% set(get(get(h(2),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
% h=area(fig_e, xx', max([res_lstp_2nd_switch_mean-res_lstp_2nd_switch_sem 2.*res_lstp_2nd_switch_sem],0));
% set(h(2),'FaceColor', black2,'EdgeColor', [1 1 1]);
% set(h(1),'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
% set(get(get(h(1),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend
% set(get(get(h(2),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'); % Exclude line from legend


plot(fig_d, xx, res_lstp_1st_switch_mean, 'Color', blue1, 'LineWidth', lineWidth);
plot(fig_d, xx, res_lstp_2nd_switch_mean, 'Color', black1, 'LineWidth', lineWidth);

if(0)
    plot(fig_d, xx, res_lstp_1st_switch_fixedU, 'Color', [0.5 0.5 1], 'LineWidth', lineWidth);
    plot(fig_d, xx, res_lstp_2nd_switch_fixedU, 'Color', [0 0 1], 'LineWidth', lineWidth);
end

set(fig_d,'layer','top')
xlabel(fig_d, 'time (s)', 'FontSize', fontsize);
ylabel(fig_d, 'p_{discrimination}', 'FontSize', fontsize+1);
legD = legend(fig_d, 'Learning', 'Relearning', 'Location', 'NorthEast');
    posD = get(legD, 'Position');
    posD(2) = posD(2)-0.04;
    set(legD, 'Position', posD);
    legend(fig_d, 'boxoff');
        
box(fig_d, 'off');
%xlim(fig_e, [-100 stime2]);
ylim(fig_d, [0.48 1.02]);
set(fig_d, 'YTick', [0.5 0.75 1])
%legE = legend(fig_e, '');
%legE = legend(fig_e, 'Learning', 'Relearning', 'Location', 'East');
%legE = legend(fig_e, 'Learning', 'Relearning', 'Learning \DeltaU=0', 'Relearning \DeltaU=0', 'Location', 'East');
%legend(fig_e, 'boxoff');
%posE = get(legE, 'Position');
%posE(1) = posE(1)+0.12;
%set(legE, 'Position', posE);
xlim(fig_d, [0 xx(end)-19.5])
set(fig_d, 'TickDir', 'Out');



%% E - Time to learn
fig_e = subplot(nrows, ncols, fig_ind(fig_subind==E), 'replace');
    pe = get(fig_e, 'Position');
    pe(1) = pe(1)+0.09;
    pe(2) = pe(2)+0.03;
    pe(3) = pe(3)-0.18;
    set(fig_e, 'Position', pe);

    %Calculate time-to-learn
    %Learning
    time_learning = zeros(size(res_lstp_1st_switch,1),1);
    for i=1:size(res_lstp_1st_switch,1)
        aux=find(res_lstp_1st_switch(i,:)>0.99);
        time_learning(i) = aux(1)./1e2;
    end
    %Relearning
    time_relearning = zeros(size(res_lstp_2nd_switch,1),1);
    for i=1:size(res_lstp_2nd_switch,1)
        aux=find(res_lstp_2nd_switch(i,:)>0.99);
        time_relearning(i) = aux(1)./1e2;
    end
    bar(fig_e, [0.5], [mean(time_learning)], 0.5, 'FaceColor', blue2, 'EdgeColor', blue1, 'LineWidth', lineWidth-0.5);
    hold on
    bar(fig_e, [1.5], [mean(time_relearning)], 0.5, 'FaceColor', black2, 'EdgeColor', black1, 'LineWidth', lineWidth-0.5);
    errorbar(fig_e, [0.5 1.5], [mean(time_learning) mean(time_relearning)], [std(time_learning)/sqrt(size(time_learning,1)) std(time_relearning)/sqrt(size(time_learning,2))], '.k', 'LineWidth', lineWidth-0.5);
    
    %Calc savings score
    savings = zeros(size(res_lstp_2nd_switch,1),1);
    for i=1:size(res_lstp_2nd_switch,1)
        savings(i) = (time_learning(i)-time_relearning(i))/time_learning(i)*100;
    end
    
    box off;
    %text(0.8, 6, ['savings score = ' num2str(mean(savings),2) ' \pm ' num2str(std(savings)/sqrt(size(savings,1)),2) '%'], 'FontSize', fontsize-1);
    xlim(fig_e, [0 2]);
    ylim(fig_e, [0 mean(time_learning)+7]);
    yLabE = ylabel(fig_e, 'time to learn (s)', 'FontSize', fontsize);
    pyLabE = get(yLabE, 'Position');
    set(yLabE, 'Position', [pyLabE(1)+0.1 pyLabE(2) pyLabE(3)]);
    set(fig_e, 'XTick', [0.5 1.5]);
    set(fig_e, 'XTickLabel', {'learn' 'relearn'});
    set(fig_e, 'TickDir', 'Out'); 
    
    
set(gcf, 'Color', 'w');

%if(save_on) export_fig(gcf, [save_path save_filename '.pdf'], '-painters','-r864', '-q101'); end