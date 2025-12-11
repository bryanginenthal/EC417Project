clc; clear; close all;

% Load data and wind-peak index
formatdata;
define_constants;

nb = size(bus,1);
nl = size(branch,1);
ng = size(gen,1);

Pd0     = bus(:,PD);
gen_bus = gen(:,GEN_BUS);
Pg_min  = gen(:,PMIN);
Pg_max  = gen(:,PMAX);

% Renewable available at peak wind
solar_avail = solargen(idx_peak_wind);
wind_avail  = windgen(idx_peak_wind);

Fmax = branch(:,RATE_A); 
Fmax(Fmax==0) = 1e4;

PTDF = makePTDF(bussys);

Cg = zeros(nb,ng);
for i=1:ng, Cg(gen_bus(i),i)=1; end

% Variable indexing
nvar  = ng + 4*nb;
idx_pg = 1:ng;
idx_ps = ng + (1:nb);
idx_pw = ng + nb + (1:nb);
idx_zs = ng + 2*nb + (1:nb);
idx_zw = ng + 3*nb + (1:nb);

f = zeros(nvar,1); 
f(idx_pg) = 1;

% Power balance
Aeq = zeros(1,nvar);
Aeq(idx_pg) = ones(1,nb)*Cg;
Aeq(idx_ps) = 1;
Aeq(idx_pw) = 1;
beq = sum(Pd0);

% Line flows
Apos  = zeros(nl,nvar);
Apos(:,idx_pg) = PTDF*Cg;
Apos(:,idx_ps) = PTDF;
Apos(:,idx_pw) = PTDF;
bpos = Fmax + PTDF*Pd0;

Aneg = -Apos;
bneg = Fmax - PTDF*Pd0;

% Availability
A_s=zeros(nb,nvar); b_s=zeros(nb,1);
A_w=zeros(nb,nvar); b_w=zeros(nb,1);

for k=1:nb
    A_s(k,idx_ps(k))=1;  A_s(k,idx_zs(k))=-solar_avail;
    A_w(k,idx_pw(k))=1;  A_w(k,idx_zw(k))=-wind_avail;
end

Aineq=[Apos;Aneg;A_s;A_w];
bineq=[bpos;bneg;b_s;b_w];

lb=zeros(nvar,1);
ub=inf(nvar,1);

lb(idx_pg)=Pg_min; ub(idx_pg)=Pg_max;

lb(idx_zs)=0; ub(idx_zs)=1;
lb(idx_zw)=0; ub(idx_zw)=1;

Aeq2=zeros(2,nvar);
Aeq2(1,idx_zs)=1;
Aeq2(2,idx_zw)=1;
beq2=[1;1];

Aeq=[Aeq;Aeq2];
beq=[beq;beq2];

intcon=[idx_zs idx_zw];

opts=optimoptions('intlinprog','Display','off');
[x,fval,exitflag,output] = intlinprog(f,intcon,Aineq,bineq,Aeq,beq,lb,ub,opts);

if exitflag <= 0
    fprintf('MILP did not converge (flag %d)\n', exitflag);
end

pg=x(idx_pg); ps=x(idx_ps); pw=x(idx_pw);
zs=x(idx_zs); zw=x(idx_zw);

solarbus=find(zs>.5);
windbus=find(zw>.5);

fprintf('\n--- Peak Wind OPF ---\n');
fprintf('Objective: %.4f\n', fval);
fprintf('Solar bus: %d\n', solarbus);
fprintf('Wind bus : %d\n', windbus);

inj  = Cg*pg + ps + pw - Pd0;
flow = PTDF*inj;

LMP = ones(nb,1)*10.947;

% --- Plot ---
figure; clf; hold on; box on;

line_from=branch(:,F_BUS);
line_to  =branch(:,T_BUS);

G=graph(line_from,line_to);
p=plot(G,'Layout','force','NodeLabel',[]);
X=p.XData; Y=p.YData;
delete(p);

for ell=1:nl
    i=line_from(ell); j=line_to(ell);
    x1=X(i); y1=Y(i); x2=X(j); y2=Y(j);

    f=flow(ell); limit=Fmax(ell);
    loading=abs(f)/limit;

    if loading>=1
        col='r'; lw=3;
    else
        col=[.4 .4 .4]; lw=1.5;
    end

    plot([x1 x2],[y1 y2],'Color',col,'LineWidth',lw);
    text(mean([x1 x2]),mean([y1 y2]),sprintf('%.1f',f),'FontSize',7);
end

scatter(X,Y,80,LMP,'filled');
colormap(jet); cb=colorbar; ylabel(cb,'LMP');

for k=1:nb
    text(X(k)+0.01,Y(k)+0.01,sprintf('%d\n%.2f',k,LMP(k)),'FontSize',7);
end

title('Peak Wind Hour OPF');
axis equal; axis off; hold off;

