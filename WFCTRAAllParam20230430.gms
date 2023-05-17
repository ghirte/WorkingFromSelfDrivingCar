$ontext
        Working From Self-Driving Cars
        Georg Hirte and Renée Laes
        2022-09-28

        Version for submission to Transport Policy

        Calculated per average worker per workday
        
        1) Benchmark: calibration for average national parameters
            (assumpotion: utility function is the same for
            each parameters combination)

            
        2) Parameter variation
            a) run for each x-parameter variation from min to max values
            (100 steps)
            b) print in files FromCar_x_Results
            
        3) Monte-Carlo simulation
         
$offtext



* ###################################################################################
* ################################################################################### 
* ###################################################################################
*               Begin DECLARATIONS of SETS and PARAMETERS
* ###################################################################################

* ====================================================================================
*====================================================================================
* Sets and Domains
* -----------------------------------------------------------------------------------

sets
i  type of contract  /a ,  b/
* a = from office, b = from car
ja(i) /a "from office"/
jb(i) /b "from car"/
l  replications of MC simulation     /1*10000/
* dynamic set used for benchmark simulation (sub(l) is used for parameter variation;
* = number of iparam without ibench and imc
subl(l) subset for parameter variation /1*10/
subl2(l) subset for MC test /1*10000/
ll(l) replications for specific simulation run 
ctr country   / DE, US /
* a = only office work, b = office + in-vehicle work
* subset - needed to refefrto specific indices in equations $-condition
subc(ctr) subitem of countries
cde(ctr) /DE /
cus(ctr) /US/
 iparam parameters  / itauw, itauf, itauc, itaus, itaud, itaup, itauq, irho,
                        igkm, igs, igd, ixbar, isp, idA, idB, irmonth,
                        isqmemployee, ibeta, ib, iy,
                        iavgcommutetime, itimeideal, iaeff, iwage, ibench, imc /
;
alias (i,j);
alias (ctr,cty);
alias(l,r1,idr);

* ====================================================================================
* ====================================================================================
* Choice of Simulation Type
* ------------------------------------------------------------------------------------
scalars
* if run_benchmark than 768: = 25
 run_benchmark     model type benchmark
* if run_parameters than 768: = <25; set to 1 if only parameter variation
 run_parameters    variation of parameters                               /1/
* if mc than 768: = 26; set to 1 if only full MC simulation
 run_mc            monte carlo simulation                                /0/
* choose country (not necessary, we always use both countries)
 run_country       MC simulation for country DE if 1, US if 0            /1/
 nloop1
 nloop2
;
* benchmark simulation (benchmark parameters)
if(run_parameters = 1, 
    nloop1 = 0;
    nloop2 = 25
elseif run_mc = 1,
    nloop1 = 25;
    nloop2 = 27;
else
    run_benchmark = 1;
    nloop1 = 24;
    nloop2 = 26;
);
* assign substrings depending on choice of simulation run
subc("DE") = yes; subc("US") = yes;

if(run_mc ne 1,
    ll(l)$(run_parameters = 1) = subl(l);
    ll("1")$(run_benchmark = 1) = yes;
else
*    ll(l) = yes;
    ll(l) = subl2(l);
);


* =====================================================================================
* =====================================================================================
* Basic Parameters
* ------------------------------------------------------------------------------------
scalar
 varphi       share of time in AV spent for leisure  /0.5/
;
parameters
 NLab(ctr)         labor demand in 1000-persons                 /DE 15000, US 70000/
 E_L0(ctr)         time endowment workday without working [h]
            /DE 8, US 8/
*++++++ benchmark parameter utility function
 y0(ctr)           option value for contract B                  /DE 0, US 0/
 aeff0(ctr)        effort parameter of effort function          /DE 2, US 2/
 delta0(ctr)       parameter utility commodity (=MUI)           /DE 1, US 1/
 delta1(ctr)       parameter utility leisure                    /DE 3.36, US 3.36/
 delta2(ctr)       parameter utility xbar (0.05)                /DE 2, US 2/
 dd2(ctr)          parameter utility xbar (ideal commute time: 16 min)  
 phi2(i,ctr)       parameter utility xbar when AV is available
 delta3(ctr)       parameter utility mobile working hours  (5)
* /DE 16, US 16/
 dd3(ctr)          parameter utility mobile working hours       
 delta1x(ctr)
 adistrib0(ctr)    parameter distribution  /DE 45, US 60/
*/DE 100, US 100/
 zA0(ctr), zA(ctr) consumption in A
 Ubar0(ctr)        fallback position: reservation utility
 indVA0(ctr)       ind utility in A
 Ubar(ctr)         fallback position: reservation utility
 indVA(ctr)        ind utility in A
 peff(ctr)         price of effort 
*++++++ time allocation
 Hwork0(ctr)       working hours per day                        /DE 8, US 8/
 Workdays0(ctr)    workdays per month                           /DE 18.2, US 19.58/
 wage0(ctr)        market wage per day (25.87 DE 25.72 US per hours)
                                                                /DE 206.96, US 205.76/
 Hworkmonth0(ctr)  working hour per month
 timeideal0(ctr)   ideal commute time in min                    /DE 16, US 16/
 avgcommutetime0(ctr) average two-way commute time in min       /DE 44, US 54.2/
 txideal0(ctr)     ideal commute time as decimals
*++++++ government
 tauw0(ctr)        wage tax average income tax rate)            /DE 0.4836, US 0.2015/
 tauf0(ctr)        add. tax per VKT private car use (EV)        /DE 0.0, US 0.0/
 tauc0(ctr)        congestion toll                              /DE 0.0, US 0.0/
 taus0(ctr)        additional tax per VKT for firm's car use    /DE 0.0, US 0.0/
 taud0(ctr)        additional tax per hour for firm's car use   /DE 0.0, US 0.0/
 taup0(ctr)        parking fee per hour with AV outside parking fee area  /DE 0.0, US 0.0/
* assume it is mainly provided by employer and private user pay only a fraction of the costs
 tauq0(ctr)        parking fee per hour at work - private car    /DE 0.5, US 0.5/
 rho0(ctr)         imputed value of the fringe benefit (private use of firm's car 30ct * 17%) /DE 0.2, US 0.2/
* travel and traffic flow
 gkm0(ctr)         net private variab monet travel costs VKT (EV BMWi3)  
 gs0(ctr)          net commercial variab monet travel costs VKT (AV) 
 gd0(ctr)          net commercial fixed travel costs per hour VKT (AV)    
 xbar0(ctr)        two-way commuting distance with office work [km]
        /DE 21.0938, US 30.687/
 avgtime0(ctr)     average commuting time  [54min -> h]         /DE 0.7285, US 0.708551/
 sp0(ctr)          share of not parking during WFC  /DE 0.5, US 0.5/
 gm0(ctr)          gross private variab travel cost VKT EV [EUR:km]  /DE 0.756519, US 0.7437832/
 gx0(ctr)          gross commercial travel costs VKT AV /DE 0.04668847, US 0.04723668/
 gh0(ctr)          gross commercial fixed travel cost per hour AV  (DE 2.56) /DE 4.005394, US 7.687163/
 tkm0(ctr)         time per km
 tb0(ctr)          parameter to determine time per km in Rietveld 1999
* ++++++ firms
 avgprod0(ctr)     avg. product per office work hour /DE 1, US 1/
 Ldemand(ctr)      labor demand
 dA0(ctr)          delay costs with contract A per day  [EUR]            /DE 1, US 1/
 dB0(ctr)          delay costs with contract B per day  [EUR]            /DE 15, US 30/
 rmonth0(ctr)      avg. rent of sqm (median for US) office space per month FFM incl. non-rental costs [EUR:sqm]
                /DE 21.391, US 91.26/
* rent per sqm /18.89/ +  other costs /3.84/ on average
 sqmemployee0(ctr) sqm of office space per employee  {sqm]               /DE 26.83, US 18.557/
 beta0(ctr)        productivity parameter per hour of in-vehicle work    /DE 0, US 0/
 b0(ctr)           fee share for private use of firm's AV on route-to-work     /DE 0, US 0/
 renth0(ctr)       avg. rent of an employee's office space per hour
* avgproduct        average productivity (value added function) per worker
 costA0(ctr)       costs per day per office worker
* help variables
 helpVTT(ctr), helpVOT(ctr)
 E_L(ctr), y(ctr), aeff(ctr),  adistrib(ctr), Hwork(ctr),  Workdays(ctr)  
 wage(ctr), Hworkmonth(ctr), timeideal(ctr), avgcommutetime(ctr), txideal(ctr), 
 tauw(ctr), tauf(ctr), tauc(ctr), taus(ctr), taud(ctr), taup(ctr), tauq(ctr),
 rho(ctr), gkm(ctr), gs(ctr), gd(ctr), xbar(ctr), avgtime(ctr)
 sp(ctr), gm(ctr), gx(ctr), gh(ctr), avgprod(ctr),
 dA(ctr), dB(ctr), rmonth(ctr), sqmemployee(ctr), beta(ctr),
 b(ctr),
 renth(ctr), costA(ctr), tb(ctr)
 ;
 b(ctr) = 0;
* =====================================================================================
* =====================================================================================
* Assignments -- calibration of utility parameters
* NOTE: We calculate time as hours with decimals instead of minutes
* -------------------------------------------------------------------------------------

* ideal commute time in hours = 16/60; 16 min as ideal commuting time in hours
txideal0(subc) = timeideal0(subc)/60;
* avgtime in hours
avgtime0(subc)  = avgcommutetime0(subc)/60;
* time per km in hours
tkm0(subc)     = max(avgtime0(subc)/xbar0(subc),0.00001);
******tkm0(subc)     = avgtime0(subc)/xbar0(subc);
* avg time per km in minutes (Rietveld et al. approach) - calculate parameter tb0 (b of Rietveld et al (1998)) for one-way trips.
tb0(subc)      = (avgcommutetime0(subc)/2 - 8 - 0.6 * max(xbar0(subc)/2-25,0))
                        /(min(xbar0(subc)/2,25));
display tb0, tkm0, xbar0, avgcommutetime0, avgtime0, txideal0;

* -+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* calibration for basic scenario without mobile work
* -------------------------------------------------------------------
* delta3 is 1/2*(8.4%+2%) of the wage (Maestas et al, 2018). WTP for working alone
delta3(subc) = 0.052*wage0(subc);
dd3(subc)    = 1;
*2/3 * Hwork0(subc);
* VTT of commuting is 0.5 (Small, 2012, Ecotra); it is about 110% higher
* than VOT of non-commuting and non-business purpose Wardman et al., TRA 2016)
* calibrate delta1 such that u_ell = VOT, we use 1. delta1 = helpVOT* (E_L - tkm*txideal)
helpVTT(subc) = 0.5*(wage0(subc)/Hwork0(subc));
helpVOT(subc) = helpVTT(subc)/1.1; 
delta1(subc) = helpVOT(subc)* (E_L0(subc) - tkm0(subc)*xbar0(subc));
* calculate parameters from VTT = VOT - u_{tx} + tau_w \rho und ideal commute time is the outcome of utlity maximization without restrictions max u(E-tx)+u(tx) -> -u_{ell}+u_{tx} = 0; 
* delta2 = VTT/(2t) . calculate the difference between VTT und VOT.
delta2(subc) = helpVTT(subc)/(2*(tkm0(subc)*xbar0(subc)-txideal0(subc)));
* dd2 = VOT + (txideal)/(txbar-taxideal)*VTT
dd2(subc) = helpVOT(subc) + (txideal0(subc))/(tkm0(subc)*xbar0(subc)-txideal0(subc)) * helpVTT(subc);
phi2("b",subc) = varphi - (1-varphi)*(helpVOT(subc)-dd2(subc))/(2*delta2(subc)*tkm0(subc)*xbar0(subc));
phi2("a",subc) = 1;
* price of effort equalizes positive utility and negative opportunity costs of additional WFC (vo)
* at 2/3 of workhours.
peff(subc) = (delta3(subc) * log(1+2/3*Hwork0(subc))) / (power(2/3*Hwork0(subc),2)); 
display helpVTT, helpVOT, E_L0, Hwork0, txideal0, delta2, dd2, delta1, phi2, peff;


* ------------------------------------------------------------------------------------
* end calibration
* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


* assign initial parameters values for simulation

E_L(ctr)=E_L0(ctr); y(ctr)=y0(ctr); aeff(ctr)=aeff0(ctr);
adistrib(ctr)=adistrib0(ctr); Hwork(ctr)=Hwork0(ctr);
Workdays(ctr) = Workdays0(ctr); wage(ctr) = wage0(ctr); Hworkmonth(ctr)=Workdays0(ctr)*Hwork0(ctr);
timeideal(ctr)=timeideal0(ctr); avgcommutetime(ctr) = avgcommutetime0(ctr);
txideal(ctr) = txideal0(ctr); tauw(ctr) = tauw0(ctr);  tauf(ctr) = tauf0(ctr);
tauc(ctr) = tauc0(ctr); taus(ctr) = taus0(ctr); taud(ctr)=taud0(ctr);
taup(ctr) = taup0(ctr);  tauq(ctr) = tauq0(ctr);
 rho(ctr) =  rho0(ctr); gkm0(ctr) = gm0(ctr)/(1+tauf0(ctr));
 gkm(ctr)= gkm0(ctr);
 gs0(ctr)=gx0(ctr)/(1+taus0(ctr)); gs(ctr) = gs0(ctr);
 gd0(ctr)=gh0(ctr)/(1+taud0(ctr)); gd(ctr) = gd0(ctr);
 xbar(ctr) =  xbar0(ctr); avgtime(ctr) = avgtime0(ctr);
 sp(ctr) =  sp0(ctr); gm(ctr) = gm0(ctr); gx(ctr)=gx0(ctr); gh(ctr) =  gh0(ctr);
 tb(ctr) =  tb0(ctr); avgprod(ctr) = avgprod0(ctr);
 dA(ctr) = dA0(ctr); dB(ctr)=dB0(ctr);
 rmonth(ctr) = rmonth0(ctr);
 sqmemployee(ctr) = sqmemployee0(ctr); beta(ctr) = beta0(ctr);
  renth(ctr) = rmonth(ctr)*sqmemployee(ctr) / Hworkmonth(ctr);
  renth0(ctr) = renth(ctr);
*b(ctr) =  b0(ctr);
 

 zA0(subc) = (1-tauw(subc))*wage(subc)
                    - ((gm(subc) + tauc(subc)*tkm0(subc))*xbar(subc))
                        - tauq(subc)*Hwork(subc) ;
 Ubar0(subc) = delta0(subc)*zA0(subc)
            + delta1(subc)*log(E_L(subc)-tkm0(subc)*xbar(subc))
                - delta2(subc)*power(tkm0(subc)*xbar(subc),2)
                    + dd2(subc)*(tkm0(subc)*xbar(subc));
Ubar(subc) = Ubar0(subc);

* ====================================================================================
* =====================================================================================
* temporary and varied parameters
* ------------------------------------------------------------------------------------
* define parameters for simulation
parameters 
mm_gkm(l), mm_gs(l), mm_gd(l), mm_xbar(l), mm_paramb(l), mm_sqmemployee(l),
mm_rmonth(l), mm_avgproduct(l), mm_wage(l),m_gkm(l), m_gs(l), m_gd(l), m_xbar(l), m_paramb(l), m_sqmemployee(l),
m_rmonth(l), m_avgproduct(l), m_wage(l),
zp_adistrib(ctr,l)        distribution parameter preference
zp_beta(ctr,l)            productivity per hour of WFC (avg daily wage +- 20 percent)
zp_gm(ctr,l)              gross monetary avg. travel costs employee [EUR:km]
zp_gx(ctr,l)              gross monetary variable travel costs firms [EUR:km]
zp_gh(ctr,l)              gross fixed variable travel costs firms per hour [EUR:h]
zp_gkm(ctr,l)
zp_gs(ctr,l)
zp_gd(ctr,l)
zp_xbar(ctr,l)            two-way commuting distance with office work (in km)
zp_avgcommutetime(ctr,l)  average commute time in min
zp_avgtime(ctr,l)         average commuting time [min -> h]
zp_sp(ctr,l)              share of riding (not parking) during in-vehicle work
zp_sqmemployee(ctr,l)     office sqm per employee
zp_dA(ctr,l)              delay costs with contract A per day  [EUR]
zp_dB(ctr,l)              delay costs with contract B per day  [EUR]
zp_costA(ctr,l)           benchmark costs of firms
zp_rmonth(ctr,l)          avg. rent of sqm office space per month FFM incl. non-rental costs [EUR:sqm]
zp_renth(ctr,l)           rents per hour of office work
zp_avgproduct(ctr,l)      average productivity (value added function) per worker
zp_timeideal(ctr,l)       ideal commute time [min]
zp_txideal(ctr,l)         ideal commute time [h]
zp_wage(ctr,l)            market wage(value added function) per worker
zp_wagestar(ctr,l)        optimal market wage (at epsilon = 0)
zp_aeff(ctr,l)            effort parameter of effort function
zp_tauw(ctr,l)            wage tax rate        
zp_tauf(ctr,l)            fuel tax rate (in addition to current rate)
zp_taus(ctr,l)            fuel tax rate commerical travel (in addition to current rate)  
zp_taud(ctr,l)            fix car cost rate
zp_taup(ctr,l)            parking fee of private car at work location
zp_tauq(ctr,l)            parking fee for AV at any location
zp_tauc(ctr,l)            congestion toll
zp_rho(ctr,l)             imputed value for private use of AV on the trip to the office
zp_tb(ctr,l)              parameter to calculate time per tkm (Rietveld 1999)
zp_b(ctr,l)               fee for private use of AV
zp_y(ctr,l)               option value for AV
zp_gamma1(ctr,l)          parameter production function 
zp_psi(ctr,l)             dimension parameter production and labor demand function
zp_tkmxbar(ctr,l)         time per km whith distance xbar
zp_phi2(i,ctr)             multiplier of utility parameter dd2 in u(tx) in B
zp_varphi                 share of VTT_B on VTT_A in empirical studies
;
zp_aeff(subc,ll)  = aeff(subc);
zp_beta(subc,ll) = beta(subc);
zp_adistrib(subc,ll) = adistrib(subc);
zp_avgproduct(subc,ll) = avgprod(subc);
zp_sp(subc,ll)    = sp(subc);
zp_dA(subc,ll)   = dA(subc);
zp_dB(subc,ll)    = dB(subc); 
zp_tauw(subc,ll) = tauw(subc); zp_tauf(subc,ll)=tauf(subc);
zp_taus(subc,ll) = taus(subc); zp_taud(subc,ll)=taud(subc);
zp_tauc(subc,ll) = tauc(subc); zp_tauq(subc,ll) = tauq(subc);
zp_taup(subc,ll) = taup(subc); zp_rho(subc,ll)= rho(subc);
zp_tb(subc,ll) = tb(subc); zp_b(subc,ll) = b(subc);
zp_y(subc,ll) = y(subc); zp_avgcommutetime(subc,ll) = avgcommutetime(subc);
zp_timeideal(subc,ll) = timeideal0(subc);
zp_txideal(subc,ll) = zp_timeideal(subc,ll)/60;
zp_xbar(subc,ll) = xbar0(subc);
zp_sqmemployee(subc,ll) = sqmemployee0(subc);
zp_wage(subc,ll) = wage0(subc);
zp_gkm(subc,ll) = gkm0(subc); zp_gm(subc,ll) = (1+zp_tauf(subc,ll))*zp_gkm(subc,ll);
zp_gs(subc,ll) = gs0(subc); zp_gx(subc,ll) = (1+ zp_taus(subc,ll))*zp_gs(subc,ll);
zp_gd(subc,ll) = gd0(subc); zp_gh(subc,ll) = (1+ zp_taud(subc,ll))*zp_gd(subc,ll);
zp_rmonth(subc,ll) = rmonth0(subc);
zp_renth(subc,ll) = zp_rmonth(subc,ll)*zp_sqmemployee(subc,ll) / Hworkmonth(subc);



* =====================================================================================
* =====================================================================================
* read in parameters for Monte-Carlo simulation (calculated in R)
* -------------------------------------------------------------------------------------
$gdxin 'MC_paramDE.gdx'
$load mm_gkm=mc_gkm mm_gs=mc_gs mm_gd=mc_gd mm_xbar=mc_xbar mm_paramb=mc_parambDE mm_sqmemployee=mc_sqmemployee mm_rmonth=mc_rmonth mm_avgproduct=mc_avgproduct mm_wage=mc_wage 
$gdxin
if( run_mc =1,
    zp_gm("DE",ll)= mm_gkm(ll);
    zp_gx("DE",ll) = mm_gs(ll); zp_gh("DE",ll)=mm_gd(ll);
    zp_xbar("DE",ll)=mm_xbar(ll);
    zp_tb("DE",ll)=mm_paramb(ll); zp_sqmemployee("DE",ll)=mm_sqmemployee(ll);
    zp_rmonth("DE",ll)=mm_rmonth(ll);
    zp_avgproduct("DE",ll)=mm_avgproduct(ll);
    zp_wage("DE",ll)=mm_wage(ll)*8;
    gm0("DE") = sum(l, mm_gkm(l))/card(l);
    gx0("DE") = sum(l, mm_gs(l))/card(l);
    gh0("DE") = sum(l, mm_gd(l))/card(l);
    xbar0("DE") = sum(l, mm_xbar(l))/card(l);
    tkm0("DE") = sum(l, mm_paramb(l))/card(l);
    sqmemployee0("DE") = sum(l, mm_sqmemployee(l))/card(l);
    rmonth0("DE") = sum(l, mm_rmonth(l))/card(l);
    avgprod0("DE") = sum(l, mm_avgproduct(l))/card(l);
    wage0("DE") = sum(l, mm_wage(l)*8)/card(l);
    wage("DE") = wage0("DE");

);

$gdxin 'MC_paramUS.gdx'
$load m_gkm=mc_gkmUS m_gs=mc_gsUS m_gd=mc_gdUS m_xbar=mc_xbarUS m_paramb=mc_parambUS m_sqmemployee=mc_sqmemployeeUS m_rmonth=mc_rmonthUS m_avgproduct=mc_avgproductUS m_wage=mc_wageUS 
$gdxin
if(run_mc =1,
    zp_gm("US",ll)= m_gkm(ll); zp_gx("US",ll)=m_gs(ll);
    zp_gh("US",ll)=m_gd(ll);
    zp_xbar("US",ll)=m_xbar(ll);
    zp_tb("US",ll)=m_paramb(ll); zp_sqmemployee("US",ll)=m_sqmemployee(ll);
    zp_rmonth("US",ll)=m_rmonth(ll); zp_avgproduct("US",ll)=m_avgproduct(ll);
    zp_wage("US",ll)=m_wage(ll)*8;
    gm0("US") = sum(l, m_gkm(l))/card(l);
    gx0("US") = sum(l, m_gs(l))/card(l);
    gh0("US") = sum(l, m_gd(l))/card(l);
    xbar0("US") = sum(l, m_xbar(l))/card(l);
    tb0("US") = sum(l, m_paramb(l))/card(l);
    sqmemployee0("US") = sum(l, m_sqmemployee(l))/card(l);
    rmonth0("US") = sum(l, m_rmonth(l))/card(l);
    avgprod0("US") = sum(l, m_avgproduct(l))/card(l);
    wage0("US") = sum(l, m_wage(l)*8)/card(l);
    wage("US") = wage0("US");
    
);


* Parameter variation in addition to random parameters in .gdx file
$funclibin stolib stodclib
Functions myduniform /stolib.duniform /
          mydtriangular / stolib.dtriangular /;
scalar ij;
if(run_mc = 1,
    loop(l(ll), zp_beta(subc,ll) =mydtriangular(-0.25,0,+0.25)* zp_wage(subc,ll));
*/Hwork0(subc));
    loop(l(ll), zp_aeff(subc,ll) = mydtriangular(1,2,3));
    loop(l(ll), zp_sp(subc,ll) = mydtriangular(0,0.5,1));
    loop(l(ll), zp_b(subc,ll) = 0);
* overhead costs for office work (only related to the office) without rents = 5% of the daily wage, ;
    loop(l(ll), zp_dA(subc,ll) = mydtriangular(10,10,10));
    loop(l(ll), zp_dB(subc,ll) = mydtriangular(5,10,15));
*    loop(l(ll), zp_y(subc,ll) = mydtriangular(0,0,0));
    zp_y(subc,ll) = 0;
);
    


* calculate parameters that are linked to randomly chosen parameters
* e.g. net travel costs per km
zp_gkm(subc,ll) = zp_gm(subc,ll)/(1+zp_tauf(subc,ll));
zp_gs(subc,ll) = zp_gx(subc,ll)/(1+zp_taus(subc,ll));
zp_gd(subc,ll) = zp_gh(subc,ll)/(1+zp_taud(subc,ll));
* monthly working hours
Hworkmonth(subc)= Hwork(subc) * workdays(subc);
* office rent per sqm and worker
zp_renth(subc,ll) = zp_rmonth(subc,ll)*zp_sqmemployee(subc,ll) / Hworkmonth(subc);
* avg. travel time in decimals calculated with Rietveld equation for one-way trip
zp_avgtime(subc,ll) = (8 +
        zp_tb(subc,ll)*min(zp_xbar(subc,ll)/2,25)
                + 0.6 * max(zp_xbar(subc,ll)/2-25,0))/60;
* avg. commuting VKT  (two-way)              
zp_tkmxbar(subc,ll) = ( zp_avgtime(subc,ll)*2 )/zp_xbar(subc,ll);
* total office costs with contract A
zp_costA(subc,ll) = Hwork(subc)*zp_renth(subc,ll) + zp_dA(subc,ll);



* ###################################################################################
*               END DECLARATIONS of SETS and PARAMETERS
* ###################################################################################
* ###################################################################################
* ################################################################################### 





* ###################################################################################
* ################################################################################### 
* ###################################################################################
*             Begin DECLARATION of EQUATIONS and MODEL
* ###################################################################################

* ====================================================================================
* ====================================================================================
* Endogenous variables
* ------------------------------------------------------------------------------------
variables
indV(ctr)                 indirect utility of contract B
indVx(ctr)                indirect utility without wage and preference parameter, contract B
u_vo(ctr), u_vovo(ctr)    1st and 2nd deriviative of utility part of vo
u_ell(ctr), u_ellell(ctr) 1st and 2nd deriviative of utility part of ell
u_tx(ctr), u_txtx(ctr)    1st and 2nd deriviative of utility part of x
epsilon(ctr)              preference for mobile work of marginal employee
elastxb(ctr)              elasticity of x wrt pb
xB(ctr)                   commuting xbar under contract B
vo(ctr)                   WFC hours substituting WFO
;

positive variables
effort(ctr)               effort to compensate for absence form office
income(ctr)               monetary income
z(ctr)                    consumption
ell(ctr)                  leisure (residual)
hvc(ctr)                  WFC hours instead of commuting(help variable)
v(ctr)                    WFC hours
alpha(ctr)                share of households chosing contract B (mobile work)
tkm (ctr)                 travel time per km
wBstar(ctr)               average net wages
xvo(ctr)                  travel distance while WFO that substitutes WFO
wB(ctr)                   reservation wage of the marginal employee 
costB(ctr)                non-wage costs per day per mobile worker
my1(ctr), myc(ctr), myv(ctr), myx(ctr)  rationing multipliers
my4(ctr), my3(ctr)        rationing parameter for pb
pb(ctr)                   private-use price of firm's AV
;


* ===================================================================================
* ===================================================================================
* Equations
* -----------------------------------------------------------------------------------
equations
Eq_effort(ctr), Eq_z_demand(ctr), Eq_ell(ctr), Eq_peff(ctr)
Eq_xB(ctr), Eq_hvc(ctr), Eq_v(ctr), Eq_vo(ctr), Eq_my1(ctr), Eq_myc(ctr), Eq_myv(ctr),
Eq_myx(ctr), Eq_alpha(ctr), Eq_xvo(ctr)
Eq_indV(ctr), Eq_indVx(ctr), 
Eq_uvo(ctr), Eq_uvovo(ctr), Eq_uell(ctr), Eq_uellell(ctr), Eq_utx(ctr), Eq_utxtx(ctr),
Eq_pb(ctr), Eq_wB(ctr),       
Eq_wBstar(ctr)        wage of marginal individual
Eq_epsilon(ctr)       preference paramter of marginal individual
Eq_costB(ctr), Eq_tkm(ctr), Eq_elastxb(ctr), Eq_my4(ctr), Eq_my3(ctr)
;

* =====================================================================================
* =====================================================================================
* Definition of equations
* -------------------------------------------------------------------------------------
* utility components (only relevant if WFC substitutes WFO)
* just to remember: subc(ctr) denotes the country
Eq_effort(ctr)$subc(ctr)..      effort(ctr) =e= vo(ctr)**aeff(ctr);
* marginal utilities
Eq_uvo(ctr)$subc(ctr)..         (vo(ctr) + dd3(ctr))*u_vo(ctr) =e= delta3(ctr);
Eq_uvovo(ctr)$subc(ctr)..       (vo(ctr) + dd3(ctr))*u_vovo(ctr) =e= - u_vo(ctr);
Eq_uell(ctr)$subc(ctr)..        u_ell(ctr)*ell(ctr) =e= delta1(ctr);
Eq_uellell(ctr)$subc(ctr)..     u_ellell(ctr)*ell(ctr) =e= -u_ell(ctr);
Eq_utxtx(ctr)$subc(ctr)..       u_txtx(ctr) =e= -2*delta2(ctr);

* FOCs
Eq_vo(ctr)$subc(ctr)..          u_vo(ctr) =e= delta0(ctr)*peff(ctr) * aeff(ctr)*(vo(ctr)**(aeff(ctr)-1))
                                        - myv(ctr) + my1(ctr);
Eq_utx(ctr)$subc(ctr)..         u_tx(ctr) =e= dd2(ctr) * phi2("b",ctr)
                                    - 2*delta2(ctr)*tkm(ctr)*xB(ctr);
Eq_xB(ctr)$subc(ctr)..          (-pb(ctr)-u_ell(ctr)+u_tx(ctr))*tkm(ctr) =e=
                                        - myx(ctr) + myc(ctr) - my1(ctr);

* subsequent. greater or equal                                        
Eq_my1(ctr)$subc(ctr)..         Hwork(ctr) - tkm(ctr)*(xbar(ctr)-xB(ctr)) - vo(ctr) =g= 0;
Eq_myv(ctr)$subc(ctr)..         vo(ctr) =g= 0;
Eq_myx(ctr)$subc(ctr)..         xB(ctr) =g= 0;
Eq_myc(ctr)$subc(ctr)..         xbar(ctr) - xB(ctr) =g= 0;

*                                              
* quantities using optimal choices of commuting distance under contract B: xB, WFC instead of WFO vo
Eq_ell(ctr)$subc(ctr)..         ell(ctr) =e= E_L(ctr)- tkm(ctr)*xB(ctr);
Eq_hvc(ctr)$subc(ctr)..         hvc(ctr) =e= tkm(ctr)*(xbar(ctr) - xB(ctr));
Eq_v(ctr)$subc(ctr)..           v(ctr) =e= hvc(ctr) + vo(ctr);
Eq_Z_demand(ctr)$subc(ctr)..    z(ctr) =e= (1-tauw(ctr))*wB(ctr)
                                    - tauw(ctr)*rho(ctr) *xbar(ctr)
                                    - pb(ctr)* tkm(ctr) * xB(ctr) 
                                    - peff(ctr)*effort(ctr);
* indirect utility
Eq_indVx(ctr)$subc(ctr)..  indVx(ctr) =e=
                                delta1(ctr)*log( ell (ctr) )
                                + dd2(ctr)*phi2("b",ctr)*( tkm(ctr)*xB(ctr) )
                                - delta2(ctr)*power(tkm(ctr)*xB(ctr),2)
                                + delta3(ctr)*log(vo(ctr)+dd3(ctr)) + y(ctr)
                                - delta0(ctr)*peff(ctr)*effort(ctr)
                                - delta0(ctr)*tauw(ctr)*rho(ctr)*xbar(ctr)
                                - delta0(ctr)*pb(ctr) * tkm(ctr) * xB(ctr) ;
Eq_indV(ctr)$subc(ctr)..   indV(ctr) =e= indVx(ctr) + (1-tauw(ctr))*wB(ctr);
Eq_wBstar(ctr)$subc(ctr)..    wBstar(ctr) * delta0(ctr) * (1-tauw(ctr))
                                    =e= (indVA(ctr) - indVx(ctr));
* end first stage employee (intensive margin and quantitity decisions)
* time per km depending on distance traveled (Rietveld approach)
Eq_tkm(ctr)$subc(ctr)..     tkm(ctr)*(max((xbar(ctr)),0.01)) =e=
                                    ( (8 + tb(ctr) * min((xbar(ctr))/2,25)
                                     + 0.6 * max((xbar(ctr))/2-25,0))*2)/60;

* we assume that travel costs per hour occurs for each hour of WFC exceeding commuting xbar
* decision of firms on payment for private use of the firm's AV
*Eq_pb(ctr)$subc(ctr)..          pb(ctr)*tkm(ctr)*(1-my4(ctr)) =e= (beta(ctr) + renth(ctr))*tkm(ctr)
*                                - xB(ctr)*(tkm(ctr)**2) * (u_ellell(ctr)+u_txtx(ctr))
*                                - my4(ctr)*(gx(ctr)+(gh(ctr)+ tauc(ctr))*tkm(ctr))
*                                - my3(ctr)*tkm(ctr)* (u_ellell(ctr)+u_txtx(ctr))
*                                ;
Eq_pb(ctr)$subc(ctr)..          pb(ctr)*tkm(ctr) =e= (beta(ctr) + renth(ctr))*tkm(ctr)
                                - (tkm(ctr)*xB(ctr)
*                                - tkm(ctr)*my4(ctr)
                                )*tkm(ctr)*
                                (u_ellell(ctr)+u_txtx(ctr))
                                ;
*Eq_my4(ctr)$subc(ctr)..         gx(ctr)+(gh(ctr)+ tauc(ctr))*tkm(ctr)
*                                 - pb(ctr)*tkm(ctr) =g= 0;
*Eq_my3(ctr)$subc(ctr)..         pb(ctr) =g= 0;  

* costs for commuting xbar: gkm/t gives travel cost per hour. sp gives share of WFC not traveled (due to parking; not yet considered: change in velocity to reduce travel costs)
Eq_xvo(ctr)$subc(ctr)..       xvo(ctr)*tkm(ctr) =e= (1-sp(ctr))*vo(ctr);

                               
Eq_costB(ctr)$subc(ctr)..     costB(ctr) =e= ( Hwork(ctr) - v(ctr) )*renth(ctr)
                                        + vo(ctr)*( (1-sp(ctr))*( gx(ctr)
                                            + tauc(ctr)*tkm(ctr) )
                                            *(1/tkm(ctr))
                                                    + sp(ctr)*taup(ctr) + gh(ctr)
                                                    )
                                        + ( gx(ctr)+(gh(ctr)
                                            + tauc(ctr))*tkm(ctr) )*xbar(ctr)
                                            + dB(ctr) - pb(ctr)* tkm(ctr)
                                            * xB(ctr) ;


* share of mobile-work contracts
Eq_epsilon(ctr)$subc(ctr)..   epsilon(ctr) =e= delta0(ctr)*(1-tauw(ctr))
                                * ( wBstar(ctr)-wage(ctr) + costB(ctr)
                                    - costA(ctr)- v(ctr)*beta(ctr) );
Eq_wB(ctr)$subc(ctr)..         wB(ctr)*(delta0(ctr)*(1-tauw(ctr))) =e=
                                wBstar(ctr)*(delta0(ctr)*(1-tauw(ctr))) - epsilon(ctr);
Eq_alpha(ctr)$subc(ctr)..     alpha(ctr) =e=
                                max( min(1/2 - (((1-tauw(ctr))
                                        *( wBstar(ctr)-wage(ctr) + costB(ctr)
                                           - costA(ctr) - v(ctr)*beta(ctr)) )
                                /(2*adistrib(ctr))),1),0);

                                    
                             



* ===================================================================================
* ===================================================================================
* Model and Algorithm
* -----------------------------------------------------------------------------------

option mcp=path;
option iterlim=900000;

* Benchmark model for basic parameter values
model WorkfromCare_Benchmark
/
  Eq_effort.effort, Eq_uvo.u_vo, Eq_uvovo.u_vovo,
  Eq_uell.u_ell, Eq_uellell.u_ellell,
  Eq_utx.u_tx, Eq_utxtx.u_txtx,
  Eq_hvc.hvc, Eq_vo.vo, Eq_tkm.tkm,
  Eq_my1.my1, Eq_myc.myc, Eq_myv.myv, Eq_myx.myx,
*  Eq_my4.my4, Eq_my3.my3,
  Eq_xB.xB, Eq_ell.ell, Eq_Z_demand.z, Eq_v.v,
  Eq_pb.pb,
  Eq_indV.indV, Eq_indVx.indVx, 
  Eq_wB.wB, Eq_xvo.xvo, Eq_costB.costB, Eq_epsilon.epsilon,
  Eq_wBstar.wBstar, Eq_alpha.alpha

/;


* ###################################################################################
*                   END DECLARATIONS and Model
* ###################################################################################
* ###################################################################################
* ################################################################################### 


* ###################################################################################
* ###################################################################################
* ################################################################################### 
*                   BEGIN BENCHMARK SIMULATION and OUTPUT
* ###################################################################################

* ===================================================================================
* ===================================================================================
* Preparing Simulations
* ----------------------------------------------------------------------------------


scalar
 rcount, rdloop, rmax, rd, s1;
 rmax = 10000;
s1=1;

parameters
* calc_... documents results of endogenous variables
calc_indV(i,ctr,l)          indirect utility of contract i
calc_uvo(ctr,l), calc_uvovo(ctr,l)  1st and 2nd deriviative of utility part of v
calc_uell(ctr,l), calc_uellell(ctr,l)        1st and 2nd deriviative of utility part of ell
calc_utx(ctr,l), calc_utxtx(ctr,l)              1st and 2nd deriviative of utility part of x
calc_income(i,ctr,l)        monetary income
calc_z(i,ctr,l)             consumption
calc_ell(i,ctr,l)           leisure (residual)
calc_xB(i,ctr,l)            commuting xbar
calc_xco(ctr,l)             additional WFC time under contract b
calc_hvc(ctr,l)             WFC worktime during commuting
calc_v(i,ctr,l)             WFC worktime
calc_vo(i,ctr,l)            WFC worktime
calc_lambda(i,ctr,l)        MUI
calc_my(i,ctr,l)            MUT
calc_my1(ctr,l), calc_myc(ctr,l), calc_myv(ctr,l), calc_myx(ctr,l),
calc_alpha(ctr,l)           share of households chosing contract b (in-vehicle work)
calc_xv(ctr,l)              travel km per hour of WFC
calc_costA(ctr,l), calc_costB(ctr,l)  costs without wage
calc_pb(ctr,l)              price for private use of firm's AV
calc_wB(i,ctr,l)            wage and omega
calc_tkm(i,ctr,l)           time per km
calc_epsilon(ctr,l)         cut-off of exogenous preference parameter
calc_wBstar(ctr,l)          wage at epsilon = 0
calc_indVx(ctr,l)           ind utility in B without net wage income
calc_effort(i,ctr,l)        effort to compensate for loss of network
calc_Ubar(ctr,l)            outside utility
calc_elastxb(ctr,l)         elasticity of x wrt b
calc_avgWTP(ctr,l)          average WTP
*calc_my3(ctr,l),
*calc_my4(ctr,l) rationing parameter pb
* zparam_... documents basic parameters that are not varied
zparam_N(ctr), zparam_E_L(ctr), zparam_delta0(ctr), zparam_delta1(ctr), zparam_delta2(ctr), zparam_dd2(ctr),
zparam_peff(ctr), zparam_delta3(ctr), zparam_dd3(ctr), zparam_Hwork(ctr),
zparam_Workdays(ctr), zparam_avgspeed(ctr,l), 
zparam_sp(ctr,l), zparam_Hworkmonth(ctr,l)
;

* Calculations using the endogenous variables 
parameters
zFOC_xb(i,ctr,l)                 check FOC
zFOC_vo(i,ctr,l)                 check FOC
zFOC_my1(i,ctr,l)                check time constraint
zemployee_VOT(i,ctr,l)           value of time
zemployee_VTT(i,ctr,l)           value of travel time savings
zemployee_ellx(i,ctr,l)          leisure
zemployee_traveltime(i,ctr,l)    commute time
zemployee_ch_budgetm(i,ctr,l)    monetary budgetconstraint of a typical employee
zemployee_ch_budgett(i,ctr,l)    time budget constraint typical employee
zemployee_ut_exputility(ctr,l)   expected utility
zemployee_ut_MUt(i,ctr,l)        MU of time per km
zemployee_ut_MUxb(i,ctr,l)       MU of travel xbar
zemployee_ut_MUv(i,ctr,l)        MU of in-vehicle time
zemployee_ut_MUell(i,ctr,l)      MU of leisure time
zemployee_ut_MUz(i,ctr,l)        MU of consumption
zemployee_ut_MUtxb(i,ctr,l)      MU of travel time
zemployee_ut_MUtxbtxb(i,ctr,l)   2nd derivative of MU of travel time
zemployee_ut_MUellell(i,ctr,l)   2nd derivative of MU of leisure (of time)
zemployee_ut_Util_ell(i,ctr,l)   utility of leisure
zemployee_ut_Util_z(i,ctr,l)     utility of z
zemployee_ut_Util_txB(i,ctr,l)     utility of commuting time tkm*x
zemployee_ut_Util_v(i,ctr,l)     utility of in-vehicle work
zemployee_tcommute(i,ctr,l)      commute time
zemployee_travelcost(i,ctr,l)    travel costs of individual
zemployee_ggen(i,ctr,l)          generalized travel costs(two-way)
zindV_z(i,ctr,l), zindV_ell(i,ctr,l), zindV_v(i,ctr,l), zindV_txB(i,ctr,l) indirect utility components
timeconstraint(ctr,l),
zVOT(i,ctr,l), zVTT(i,ctr,l), zalpha_calc(ctr,l),

firm_production(i,ctr,l)        productivity function
firm_margprod(i,ctr,l)
firm_avg_expwage(i,ctr,l)       average expected costs per worker
firm_NetProfits(i,ctr,l)        Net profits of firms
firm_labor_demand(i,ctr,l)      labor demand
firm_labor_change(ctr,l)        change in labor demand
firm_FOC_labordemand(i,ctr,l)   Check whether FOC for labor demand is fullfilled
firm_chg_margprod(ctr,l)        change in marginal productivity by moving to non-office worker
firm_chg_travel(ctr,l)            check travel costs of firms per day and worker type
firm_marg_Deltafc(ctr,l)        non-wage fix cost difference between B and A
firm_marg_Deltamc(ctr,l)        non-wage marginal cost difference between B and A (per h of mobile work)
firm_marg_Deltac(ctr,l)         change in non-wage costs if hiring a marginal non-office worker
firm_marg_costsA(ctr,l)         marginal costs contract A per day and office worker (no wage)
firm_marg_costsB(ctr,l)         marginal costs contract B per day and non-office worker (no wage)
firm_chg_marg_costs(ctr,l)      change in costs per day if hiring a marginal mobile worker
firm_chg_marg_profit(ctr,l)     change in profits per working hour of a typical firm
elast_labdemand(ctr,l)          elasticity of labor demand
firm_avg_sqm(i,ctr,l)           avgerage sqm per worker corrected for mobile work
firm_sqm(i,ctr,l)               sum of sqm demand by firm
firm_chg_avg_sqm(i,ctr,l)       percentage change in avg. sqm per worker
firm_chg_sqm(ctr,l)             percentage change in sqm demand of firm
;

zparam_delta1(ctr) = delta1(ctr);
zparam_delta2(ctr) = delta2(ctr);
zparam_delta3(ctr) = delta3(ctr);
zparam_dd2(ctr) = dd2(ctr);
zparam_dd3(ctr) = dd3(ctr);


* declarations for output files
* set -> gives output file names
set
gdxFiles Outputnamen der gdx Files       /MobileWork_Benchmark,MobileWork_tauw,MobileWork_tauf,
                                        MobileWork_tauc,MobileWork_taus,MobileWork_taud,MobileWork_taup,
                                        MobileWork_tauq, MobileWork_rho, MobileWork_gm, MobileWork_gx,
                                        MobileWork_gh, MobileWork_xbar, MobileWork_sp, MobileWork_dA,
                                        MobileWork_dB, MobileWork_rmonth, MobileWork_sqmemployee,
                                        MobileWork_beta, MobileWork_b, MobileWork_y, MobileWork_tb,
                                        MobileWork_tideal, MobileWork_eff, MobileWork_wage,
                                        MobileWork_MonteCarlo
                                        /
s        test set          /s1*s26/
;

* assigns number to the outputfiles, used as index
parameter iford(gdxFiles)  Dummy fuer ord(iparam) /MobileWork_Benchmark=25,MobileWork_tauw=1,MobileWork_tauf=2,
            MobileWork_tauc=3,MobileWork_taus=4,MobileWork_taud=5,MobileWork_taup=6,
            MobileWork_tauq=7, MobileWork_rho=8, MobileWork_gm=9, MobileWork_gx=10, MobileWork_gh=11,
            MobileWork_xbar=12, MobileWork_sp=13, MobileWork_dA=14, MobileWork_dB=15, MobileWork_rmonth=16,
            MobileWork_sqmemployee=17, MobileWork_beta=18, MobileWork_b=19, MobileWork_y=20, MobileWork_tb=21,
            MobileWork_tideal=22, MobileWork_eff=23, MobileWork_wage=24, MobileWork_MonteCarlo=26/;
            
scalar aa/10/, bb/20/;


* =======================================================================================
* =======================================================================================
* Start outer Simulation and Outputloop
* ---------------------------------------------------------------------------------------
scalar
 namefile, test1, check1, numberloop;
acronym aparam;
parameter oiparam;

rd = card(ll);
* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*display zp_tauw, zp_gkm, calc_tkm, zp_avgcommutetime, zp_avgtime, zp_psi, zp_xbar;
* loop over parameters - basic values

loop(iparam$(ord(iparam) > nloop1 and ord(iparam) < nloop2),
test1 = ord(iparam);
display test1, run_mc, run_parameters, run_benchmark;

* assign basic parameters to the parameter names used in equations
    beta(subc) = beta0(subc);
    gm(subc) = gm0(subc);
    gx(subc) = gx0(subc);
    gh(subc) = gh0(subc);
    gs(subc) = gs0(subc);
    gkm(subc) = gkm0(subc);
    gd(subc)  = gd0(subc);
    xbar(subc) = xbar0(subc);
    avgtime(subc) = avgtime0(subc);
    sp(subc) = sp0(subc);
    sqmemployee(subc) = sqmemployee0(subc);
    dA(subc) = dA0(subc);
    dB(subc) = dB0(subc);
    rmonth(subc) = rmonth0(subc);
    renth(subc) = renth0(subc);
    wage(subc) = wage0(subc);
    aeff(subc) = aeff0(subc);
    tauw(subc) = tauw0(subc);
    tauf(subc) = tauf0(subc);
    taus(subc) = taus0(subc);
    taud(subc) = taud0(subc);
    tauc(subc) = tauc0(subc);
    taup(subc) = taup0(subc) ;
    tauq(subc) = tauq0(subc) ;
    tb(subc) = tb0(subc);
    rho(subc) = rho0(subc) ;
    avgprod(subc) = avgprod0(subc) ;
    adistrib(subc) = adistrib0(subc) ;
    y(subc) = y0(subc);
    b(subc) = b0(subc);
    timeideal(subc) = timeideal0(subc);
    txideal(subc) = timeideal(subc)/60;

   
* Assign parameter values 
* define start values for each loop
   if( ord(iparam) eq 1,
            zp_tauw(subc,ll) = 0;
   elseif( ord(iparam) eq 2), zp_tauf(subc,ll) = tauf0(subc);
                    zp_gm(subc,ll) =
                        ((1+zp_tauf(subc,ll))*zp_gkm(subc,ll)); 
   elseif( ord(iparam) eq 3), zp_tauc(subc,ll) = tauc0(subc); 
   elseif( ord(iparam) eq 4), zp_taus(subc,ll) = taus0(subc);
*(taus0(subc)-0.2);
                zp_gx(subc,ll) = ((1+zp_taus(subc,ll))*zp_gs(subc,ll)); 
   elseif( ord(iparam) eq 5), zp_taud(subc,ll) = taud0(subc);
*(taud0(subc)-0.5);
                zp_gh(subc,ll) = ((1+zp_taud(subc,ll))*zp_gd(subc,ll)); 
   elseif( ord(iparam) eq 6), zp_taup(subc,ll) = taup0(subc); 
   elseif( ord(iparam) eq 7), zp_tauq(subc,ll) = tauq0(subc); 
   elseif( ord(iparam) eq 8), zp_rho(subc,ll) = 0; 
   elseif( ord(iparam) eq 9), zp_gkm("DE",ll) = 0.279;
                        zp_gkm("US",ll) = 0.1931;
                        zp_gm(subc,ll) = ((1+zp_tauf(subc,ll))*zp_gkm(subc,ll));
   elseif( ord(iparam) eq 10), zp_gs("DE",ll) = 0.038;
                        zp_gs("US",ll) = 0.0157;
                    zp_gx(subc,ll) = ((1+zp_taus(subc,ll))*zp_gs(subc,ll)); 
   elseif( ord(iparam) eq 11), zp_gd("DE",ll) = 1.04;
                          zp_gd("US",ll) = 7.0157; 
                    zp_gh(subc,ll) = ((1+zp_taud(subc,ll))*zp_gd(subc,ll)); 
   elseif( ord(iparam) eq 12), zp_xbar(subc,ll) = 2;
            calc_tkm('a',subc,ll) = ((8 + zp_tb(subc,ll) *
                                    min(zp_xbar(subc,ll)/2,25)
                                + 0.6 * max(zp_xbar(subc,ll)/2-25,0))*2)
                                / zp_xbar(subc,ll) / 60;
            zp_avgtime(subc,ll) = calc_tkm('a',subc,ll) * zp_xbar(subc,ll);
            zp_avgcommutetime(subc,ll) = zp_avgtime(subc,ll)*60; 
   elseif( ord(iparam) eq 13), zp_sp(subc,ll) = 0; 
   elseif( ord(iparam) eq 14), zp_dA(subc,ll) = 1; 
   elseif( ord(iparam) eq 15), zp_dB(subc,ll) = 1; 
   elseif( ord(iparam) eq 16), zp_rmonth("DE",ll) = 8.44;
                          zp_rmonth("US",ll) = 22.79;
            zp_renth(subc,ll) =
                    ((zp_rmonth(subc,ll)*zp_sqmemployee(subc,ll))/ Hworkmonth(subc)); 
   elseif( ord(iparam) eq 17),  zp_sqmemployee("DE",ll) = 18.8;
                                zp_sqmemployee("US",ll) = 9.3; 
                zp_renth(subc,ll) =
                        ((zp_rmonth(subc,ll)*zp_sqmemployee(subc,ll))/ Hworkmonth(subc)); 
   elseif( ord(iparam) eq 18), zp_beta(subc,ll) = -0.25*wage0(subc); 
* -0.25   
   elseif( ord(iparam) eq 19), zp_b(subc,ll) = 0; 
   elseif( ord(iparam) eq 20), zp_y(subc,ll) = 0; 
   elseif( ord(iparam) eq 21), zp_tb(subc,ll) = 1;
*           zp_avgcommutetime(subc,ll) = 
*            zp_avgtime(subc,ll) = calc_tkm(subc,ll) * zp_xbar(subc,ll);
*            zp_avgcommutetime(subc,ll) = zp_avgtime(subc,ll)*60; 
   elseif( ord(iparam) eq 22), zp_timeideal(subc,ll) = 5; 
                                zp_txideal(subc,ll) = (zp_timeideal(subc,ll)/60); 
   elseif( ord(iparam) eq 23), zp_aeff(subc,ll) = 1; 
   elseif( ord(iparam) eq 24), zp_wage("DE",ll) = 102.32;
                        zp_wage("US",ll) = 80.56;
    );

* =======================================================================================
* =======================================================================================
* Start Inner Simulation Loop
* ---------------------------------------------------------------------------------------
* loop for the chosen parameter
    rdloop = 0;
    loop(ll$(rdloop <= rd),
        rdloop = rdloop +1;
        display ll, rdloop, rd;
* new parameter - updated in each run
        if(ord(iparam) eq 1, 
            tauw(subc)  = zp_tauw(subc,ll) + 0.075*(rdloop-1);
        elseif(ord(iparam) eq 2),
            tauf(subc)  = zp_tauf(subc,ll) + (0.1*(rdloop-1));
            gm(subc) = ((1+tauf(subc))*gkm(subc));
        elseif(ord(iparam) eq 3),
            tauc(subc)  = zp_tauc(subc,ll) + (1*(rdloop-1));
* tauc: max 5 EUR, avg. commuting travel time is 22 min., we set a maximum of 10 EUR per h, 
        elseif(ord(iparam) eq 4),
            taus(subc)  = zp_taus(subc,ll) + (0.1*(rdloop-1));
            gx(subc) = (gs(subc) * (1+taus(subc)));
        elseif(ord(iparam) eq 5),
            taud(subc)  = zp_taud(subc,ll) + (0.1*(rdloop-1)); 
            gh(subc) = ((1+taud(subc)) * gd(subc)); 
        elseif(ord(iparam) eq 6),
            taup(subc)  = zp_taup(subc,ll) + (0.5*(rdloop-1));
* parking fee per hour anyplace            
        elseif(ord(iparam) eq 7),
            tauq(subc)  = zp_tauq(subc,ll) + (0.5*(rdloop-1));
* parking fee per hour near work
        elseif(ord(iparam) eq 8),
            rho(subc)   = zp_rho(subc,ll) + (0.1*(rdloop-1));
        elseif(ord(iparam) eq 9), 
            gkm("DE") = zp_gkm("DE",ll) + 0.09*(rdloop-1);
            gkm("US") = zp_gkm("US",ll) + 0.103*(rdloop-1);
            gm(subc) = ((1+tauf(subc))*gkm(subc)); 
        elseif(ord(iparam) eq 10),
            gs("DE") = zp_gs("DE",ll) + 0.0017*(rdloop-1);
            gs("US") = zp_gs("US",ll) + 0.0055*(rdloop-1);
            gx(subc) = (gs(subc) * (1+taus(subc)));
        elseif(ord(iparam) eq 11),
            gd("DE") = zp_gd("DE",ll) + 0.6*(rdloop-1);
            gd("US") = zp_gd("US",ll) + 0.0878*(rdloop-1);
            gh(subc) = ((1+taud(subc)) * gd(subc)); 
        elseif(ord(iparam) eq 12),
            xbar("DE") = zp_xbar("DE",ll) + 10*(rdloop-1);
            xbar("US") = zp_xbar("US",ll) + 15*(rdloop-1);
        elseif(ord(iparam) eq 13),
            sp(subc) = zp_sp(subc,ll) + (0.05*(rdloop-1)); 
        elseif(ord(iparam) eq 14),
            dA(subc) = zp_dA(subc,ll) + (1*(rdloop-1));
        elseif(ord(iparam) eq 15),
            dB(subc) = zp_dB(subc,ll) + (3*(rdloop-1));
        elseif(ord(iparam) eq 16),
            rmonth("DE") = zp_rmonth("DE",ll) + 3.68*(rdloop-1);
            rmonth("US") = zp_rmonth("US",ll) + 13.7*(rdloop-1);
            renth(subc) = ((rmonth(subc)*sqmemployee(subc))
                        / Hworkmonth(subc)); 
        elseif(ord(iparam) eq 17),
            sqmemployee("DE") = zp_sqmemployee("DE",ll)
                                + 1.6*(rdloop-1);
            sqmemployee("US") = zp_sqmemployee("US",ll)
                                + 1.86*(rdloop-1); 
            renth(subc) = ((rmonth(subc)*sqmemployee(subc))
                        / Hworkmonth(subc)); 
        elseif(ord(iparam) eq 18),
            beta(subc) = zp_beta(subc,ll) + (0.05*(rdloop-1))*wage0(subc); 
*         beta(subc) = zp_beta(subc,ll) + (0.05*(rdloop-1));       
        elseif(ord(iparam) eq 19),
            b(subc) = zp_b(subc,ll) + (0.5*(rdloop-1));
*+ (0.1*(rdloop-1)); 
        elseif(ord(iparam) eq 20),
            y(subc) = zp_y(subc,ll) + (0.5*(rdloop-1)); 
        elseif(ord(iparam) eq 21),
            tb(subc) = tb(subc) + 0.2 *(rdloop-1);
*            avgtime(subc)  = tkm(subc)*xbar(subc);
*            avgcommutetime(subc) = avgtime(subc) * 60; 
        elseif(ord(iparam) eq 22),
            timeideal(subc) = zp_timeideal(subc,ll) +  (1*(rdloop-1)); 
            txideal(subc) = (timeideal(subc )/60); 
        elseif(ord(iparam) eq 23),
            aeff(subc) = zp_aeff(subc,ll);
            if (rdloop >1, aeff(subc) = zp_aeff(subc,ll) +  0.1*rdloop ); 
        elseif(ord(iparam) eq 24),
            wage("DE") = zp_wage("DE",ll) + 50 *(rdloop-1);
            wage("US") = zp_wage("US",ll) + 50*(rdloop-1);
        elseif(ord(iparam) eq 26),
            beta(subc) = zp_beta(subc,ll);
            gm(subc) = zp_gm(subc,ll);
            gx(subc) = zp_gx(subc,ll);
            gh(subc) = zp_gh(subc,ll);
            gs(subc) = zp_gs(subc,ll);
            gkm(subc) = zp_gkm(subc,ll);
            gd(subc)  = zp_gd(subc,ll);
            xbar(subc) = zp_xbar(subc,ll);
            avgtime(subc) = avgtime0(subc);
            sp(subc) = zp_sp(subc,ll);
            sqmemployee(subc) = zp_sqmemployee(subc,ll);
            dA(subc) = zp_dA(subc,ll);
            dB(subc) = zp_dB(subc,ll);
            rmonth(subc) = zp_rmonth(subc,ll);
            renth(subc) = zp_renth(subc,ll);
            wage(subc) = zp_wage(subc,ll);
            aeff(subc) = zp_aeff(subc,ll);
            avgprod(subc) = zp_avgproduct(subc,ll);
            tb(subc) = zp_tb(subc,ll) ;
            y(subc) = zp_y(subc,ll);
*            b(subc) = zp_b(subc,ll);
        );
    

* calculate parameters  
        costA(subc) = Hwork(subc) * renth(subc) + dA(subc);
        avgprod(subc) = costA(subc) + wage(subc);
* consumption in A is fixed (no endogenous variables) (two-way)
        zp_tkmxbar(subc,ll) = (( 8 + tb(subc) * min(xbar(subc)/2,25)
                         + 0.6 * max( xbar(subc)/2 - 25,0) )*2)
                         / 60 /(xbar(subc));
        zA(subc)    = (1-tauw(subc))*wage(subc)
                    - ((gm(subc) + tauc(subc)*zp_tkmxbar(subc,ll))*xbar(subc))
                        - tauq(subc)*Hwork(subc);
* outside utility is equal to indirect utility in A
        indVA(subc) = delta0(subc)*zA(subc)
            + delta1(subc)*log(E_L(subc)-zp_tkmxbar(subc,ll)*xbar(subc))
            + dd2(subc)*zp_tkmxbar(subc,ll)*xbar(subc)
            - delta2(subc)*power(zp_tkmxbar(subc,ll)*xbar(subc),2);
        
  
* Begin set parameters     
*    Begin Starting Values
        z.l(subc)= 50;
        xb.l(subc) = 1;
        hvc.l(subc) = 0.55;
        vo.l(subc) = 1.0;
        wB.l(subc) = wage(subc);
        epsilon.l(subc) = 0.1;
        wBstar.l(subc) = wage(subc)*0.8;
        alpha.l(subc) = 0;
        myx.l(subc) = 0; myv.l(subc) = 0; my1.l(subc) = 1; myc.l(subc)=0;
        pb.l(subc) = 0.5;
        tkm.l(subc) = tkm0(subc);
*        my3.l(subc) = 0.0;
*        my4.l(subc) = 0.5;
*scalar line ;


* define lower and upper bounds
*xb.lo(i) = 0;
        z.lo(ctr) = 0.0000011;
        ell.lo(ctr) = 0.01;
        ell.up(ctr) = E_L(ctr);
***        xb.lo(ctr) = 0;
***        xb.up(ctr) = xbar(ctr);
        tkm.lo(ctr) = 0.001;
        tkm.up(ctr) = 500;
        xB.lo(ctr) = 0;
        xB.up(ctr) = xbar(ctr);
*hvo.lo = 0;
        hvc.lo(ctr) = 0;
        v.lo(ctr) = 0;
***        vo.lo(ctr) = 0;
        alpha.lo(ctr) = 0;
        alpha.up(ctr) = 1;
        pb.lo(ctr) = 0;
        pb.up(ctr) = gx(ctr)/tkm.l(ctr)+gh(ctr)+ tauc(ctr);
        wB.lo(ctr) = 5;
        wBstar.lo(ctr) =5;
        u_vovo.l(ctr) = -0.02;
        epsilon.lo(ctr) = - adistrib(ctr);
        epsilon.up(ctr) = adistrib(ctr);



* run MC-calculations for model
        for(s1 = 1 to 1 by 1000,
* End starting values (Define a loop to adjust starting values
            display s1;
            solve WorkfromCare_Benchmark using mcp;
            if ((WorkfromCare_Benchmark.modelstat =1),s1=s1+100000000000 );
        );


* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Assign Simulation Results to Output Parameters
* --------------------------------------------------------------------
        calc_tkm("b",subc,ll) = tkm.l(subc);
        calc_tkm("a",subc,ll) = zp_tkmxbar(subc,ll);
        calc_indV(i,subc,ll) = indVA(subc)$(ja(i)) + indV.l(subc)$(jb(i)); 
        calc_uvo(subc,ll) = u_vo.l(subc);
        calc_uvovo(subc,ll) = u_vovo.l(subc);
        calc_uell(subc,ll) = u_ell.l(subc);
        calc_uellell(subc,ll) = u_ellell.l(subc);
        calc_utx(subc,ll) = u_tx.l(subc);
        calc_utxtx(subc,ll) = u_txtx.l(subc);
        calc_z(i,subc,ll) = zA(subc)$(ja(i)) + z.l(subc)$(jb(i));
        calc_ell(i,subc,ll) = (E_l(subc) - calc_tkm(i,subc,ll)*xbar(subc))$(ja(i))
                                + (E_l(subc) - calc_tkm(i,subc,ll)*xb.l(subc))$(jb(i));
        calc_xb(i,subc,ll) = xbar(subc)$(ja(i)) + xb.l(subc)$(jb(i));
        calc_hvc(subc,ll) = hvc.l(subc);
        calc_v(i,subc,ll) = v.l(subc)$(jb(i));
        calc_vo(i,subc,ll) = vo.l(subc)$(jb(i));
        calc_my(i,subc,ll) = (delta1(subc)/(E_L(subc)-
                            calc_tkm(i,subc,ll)*xbar(subc)))$(ja(i))
                    + (delta1(subc)/(E_L(subc) - calc_tkm(i,subc,ll)*xB.l(subc)))$(jb(i)) ;
        calc_my1(subc,ll) = my1.l(subc);
        calc_myc(subc,ll) = myc.l(subc);
        calc_myv(subc,ll) = myv.l(subc);
        calc_myx(subc,ll) = myx.l(subc);
*        calc_my3(subc,ll) = my3.l(subc);
*        calc_my4(subc,ll) = my4.l(subc);
        calc_alpha(subc,ll) = alpha.l(subc); 
        calc_xv(subc,ll)    = xvo.l(subc);
        calc_pb(subc,ll) = pb.l(subc);
        calc_wB(i,subc,ll) = wage(subc)$(ja(i))  + wB.l(subc)$(jb(i));
        calc_epsilon(subc,ll) = epsilon.l(subc);
        calc_wBstar(subc,ll)= wBstar.l(subc);
        calc_indVx(subc,ll) = indVx.l(subc);
        calc_costA(subc,ll) = costA(subc);
        calc_costB(subc,ll) = costB.l(subc);
        calc_effort(i,subc,ll) = effort.l(subc)$(jb(i));
        calc_Ubar(subc,ll) = indVA(subc);
        calc_epsilon(subc,ll) = epsilon.l(subc);
        calc_elastxb(subc,ll) = - calc_pb(subc,ll) / (
                                (calc_xB("b",subc,ll)$(calc_xB("b",subc,ll)>0)
                                 + 0.1$(calc_xB("b",subc,ll)=0))*
                                calc_tkm("b",subc,ll)
                                * (calc_uellell(subc,ll)+calc_utxtx(subc,ll)));
                                
        zparam_E_L(subc)=E_L(subc); zparam_Hwork(subc) = Hwork(subc);
        zparam_Workdays(subc) = Workdays(subc);
        zparam_Hworkmonth(subc,ll) = Hworkmonth(subc);
        zparam_peff(subc) = peff(subc);
        zp_beta(subc,ll) = beta(subc);
        zp_gm(subc,ll) = gm(subc);
        zp_gx(subc,ll) = gx(subc);
        zp_gh(subc,ll) = gh(subc);
        zp_gs(subc,ll) = gs(subc);
        zp_gkm(subc,ll) = gkm(subc);
        zp_gd(subc,ll)  = gd(subc);
        zp_xbar(subc,ll) = xbar(subc);
        zp_tb(subc,ll) = tb(subc);
        zp_avgtime(subc,ll) = avgtime(subc);
        zp_sp(subc,ll) = sp(subc);
        zp_sqmemployee(subc,ll) = sqmemployee(subc);
        zp_dA(subc,ll) = dA(subc);
        zp_dB(subc,ll) = dB(subc);
        zp_rmonth(subc,ll) = rmonth(subc);
        zp_renth(subc,ll) = renth(subc); 
        zp_aeff(subc,ll) = aeff(subc);
        zp_tauw(subc,ll) = tauw(subc);
        zp_tauf(subc,ll) = tauf(subc);
        zp_taus(subc,ll) = taus(subc);
        zp_taud(subc,ll) = taud(subc);
        zp_tauc(subc,ll) = tauc(subc);
        zp_taup(subc,ll) = taup(subc);
        zp_tauq(subc,ll) = tauq(subc);
        zp_rho(subc,ll) = rho(subc);
        zp_avgproduct(subc,ll) = avgprod(subc);
        zp_adistrib(subc,ll) = adistrib(subc);
        zp_sp(subc,ll) = sp(subc);
        zp_y(subc,ll) = y(subc);
        zp_b(subc,ll) = b(subc);
        zp_costA(subc,ll) = costA(subc);
        zp_wage(subc,ll) = wage(subc);
        zp_phi2(i,subc) = phi2(i,subc);
        zp_varphi = varphi;
        zp_timeideal(subc,ll) = timeideal(subc);
        zp_txideal(subc,ll) = txideal(subc);
        

* calculate employment
*    gamma1 * log(employ(subc,ll)) + calc_alpha(subc,ll)*zp_beta(subc,ll)*n - n* ((1-alpha)*wage + alpha*wb) -  employ(subc,ll) * (calc_alpha(subc,ll)*calc_costB(subc,ll) + (1-calc_alpha(subc,ll)*calc_costA(subc,ll));  
* employment per firm
   
    
    );
* ---------------------------------------------------------------------------------------
* end simulation loop over index ll 
* =======================================================================================
* =======================================================================================


* ==========================================================================================
* ==========================================================================================
* Final Calculations: Check Constraints and further Calculations
* ------------------------------------------------------------------------------------------


*zparam_flow0 = flow0;  zparam_avgspeed = avgspeed; zparam_f0 = f0;
*zparam_f1 = f1; zparam_f2 = f2; zparam_capacity = capacity;
    zparam_delta0(subc) = delta0(subc);
    zparam_delta1(subc) = delta1(subc);
    zparam_delta2(subc) = delta2(subc);
    zparam_dd2(subc) = dd2(subc); zparam_delta3(subc) = delta3(subc);
* marginal utilities
    zindV_z(i,subc,ll)           = zparam_delta0(subc) * calc_z(i,subc,ll);
    zindV_ell(i,subc,ll)         = zparam_delta1(subc)*log(zparam_E_L(subc)
                                    -calc_tkm(i,subc,ll)*calc_xB(i,subc,ll));
    zindV_v(i,subc,ll)           = zparam_delta3(subc)*log(calc_vo(i,subc,ll)+1);
    zindV_txB(i,subc,ll)         = - zparam_delta2(subc)*power(calc_tkm(i,subc,ll)
                                *calc_xb(i,subc,ll),2)
                                + zparam_dd2(subc)*phi2(i,subc)*(calc_tkm(i,subc,ll)
                                    *calc_xb(i,subc,ll));
    zemployee_ut_MUtxb(i,subc,ll) = - 2*zparam_delta2(subc)*calc_tkm(i,subc,ll)
                                * calc_xb(i,subc,ll) + zparam_dd2(subc)*phi2(i,subc);
    zemployee_ut_MUt(i,subc,ll)   = - (2*zparam_delta2(subc)*calc_tkm(i,subc,ll)
                                    *calc_xb(i,subc,ll)
                                    - zparam_dd2(subc)*phi2("b",subc))*calc_xb(i,subc,ll);                               
    zemployee_ut_MUxb(i,subc,ll)  = - (2*zparam_delta2(subc)*calc_tkm(i,subc,ll)
                                    *calc_xb(i,subc,ll)
                                    - zparam_dd2(subc)*phi2(i,subc))*calc_tkm(i,subc,ll); 
    zemployee_ut_MUv(i,subc,ll)   = zparam_delta3(subc)
                                    /(calc_v(i,subc,ll) - dd3(subc));
    zemployee_ut_MUell(i,subc,ll) = zparam_delta1(subc) /(zparam_E_L(subc)
                                        -calc_tkm(i,subc,ll)*calc_xb(i,subc,ll));
    zemployee_ut_MUz(i,subc,ll)   = zparam_delta0(subc);
    zemployee_VOT(i,subc,ll)      = zemployee_ut_MUell(i,subc,ll)/zparam_delta0(subc);
    zemployee_tcommute(i,subc,ll) = (calc_tkm(i,subc,ll)*zp_xbar(subc,ll))$(ja(i))
                                 + (calc_tkm(i,subc,ll)*calc_xB(i,subc,ll))$(jb(i));
*/calc_lambda(i,l);
* control of solution: budgetm and budgett must be zero 
    zemployee_travelcost("a",subc,ll) = (zp_gm(subc,ll)
                    + zp_tauc(subc,ll)*calc_tkm("a",subc,ll))
                    * calc_xB("a",subc,ll)  + zp_tauq(subc,ll)*Hwork(subc);
    zemployee_travelcost("b",subc,ll) =
                    zp_tauw(subc,ll)*zp_rho(subc,ll)*zp_xbar(subc,ll) 
                    - calc_pb(subc,ll)*calc_tkm("b",subc,ll)*calc_xB("b",subc,ll);
    calc_income(i,subc,ll) = (1-zp_tauw(subc,ll))*( calc_wB("b",subc,ll)$jb(i)
                            + calc_wB("a",subc,ll)$ja(i) )
                         - zemployee_travelcost(i,subc,ll);
    zemployee_ch_budgetm("a",subc,ll)  = (1-zp_tauw(subc,ll))
                            * calc_wB("a",subc,ll) - (zp_gm(subc,ll)
                    + zp_tauc(subc,ll)*calc_tkm("b",subc,ll))
                    * calc_xB("a",subc,ll)  - zp_tauq(subc,ll)*Hwork(subc)
                                        - calc_z("a",subc,ll);
                                        
    zemployee_ch_budgetm("b",subc,ll)  = (1-zp_tauw(subc,ll)) * calc_wB('b',subc,ll) -
                    zp_tauw(subc,ll)*zp_rho(subc,ll)* zp_xbar(subc,ll)
                        - calc_pb(subc,ll)* calc_tkm("b",subc,ll)
                        * calc_xB("b",subc,ll) - calc_z("b",subc,ll);
    zemployee_traveltime(i,subc,ll)    =  calc_tkm(i,subc,ll) * calc_xb(i,subc,ll);
    zemployee_ch_budgett(i,subc,ll)    = zparam_E_L(subc) - calc_ell(i,subc,ll)
                                            - calc_tkm(i,subc,ll) * calc_xb(i,subc,ll);
    zemployee_ut_exputility(subc,ll)   =  (1-calc_alpha(subc,ll))* calc_indV('a',subc,ll)
                                            + calc_alpha(subc,ll) * calc_indV('b',subc,ll);
    zemployee_ggen(i,subc,ll)  = calc_uell(subc,ll)*calc_tkm(i,subc,ll)
                                + zp_gm(subc,ll);

    zVOT(i,subc,ll)    = zparam_delta1(subc)/(zparam_E_L(subc)
                                - calc_tkm(i,subc,ll)*calc_xb(i,subc,ll)) / delta0(subc);
    zVTT(i,subc,ll)    = zVOT(i,subc,ll) + (2*zparam_delta2(subc)
                                *calc_tkm(i,subc,ll) *calc_xb(i,subc,ll)
                                - zparam_dd2(subc)*phi2(i,subc))/ delta0(subc);                

    zFOC_vo(i,subc,ll) = calc_uvo(subc,ll) - delta0(subc)*zparam_peff(subc)
                                * zp_aeff(subc,ll)*
                                (calc_vo(i,subc,ll)**(zp_aeff(subc,ll)-1))
                                + calc_myv(subc,ll) - calc_my1(subc,ll) ;
    zFOC_xb(i,subc,ll) = -(calc_uell(subc,ll)+calc_utx(subc,ll))
                                *calc_tkm(i,subc,ll) - delta0(subc)*zp_tauw(subc,ll)
                                *zp_rho(subc,ll)
                                + calc_myx(subc,ll) - calc_myc(subc,ll)
                                + calc_my1(subc,ll);
    zFOC_my1(i,subc,ll) = (Hwork(subc) - calc_tkm(i,subc,ll)
                                    *(zp_xbar(subc,ll)
                                -calc_xB(i,subc,ll) ) - calc_vo(i,subc,ll))/delta0(subc);
    zalpha_calc(subc,ll) = 1/2 - ( (1-zp_tauw(subc,ll))*( calc_wB('b',subc,ll)
                                        -calc_wB('a',subc,ll)
                            + calc_costB(subc,ll) - calc_costA(subc,ll) -
                                calc_v('b',subc,ll)*zp_beta(subc,ll)) )
                            /(2*zp_adistrib(subc,ll));


* calculations marginal decision of firms: hiring a mobile or office worker
    firm_marg_Deltafc(subc,ll)       = ( ((1+zp_taud(subc,ll))*zp_gd(subc,ll)
                                 + zp_tauc(subc,ll))*calc_tkm('b',subc,ll))
                                 * zp_xbar(subc,ll)
                                 + zp_dB(subc,ll) - zp_dA(subc,ll);
    firm_marg_Deltamc(subc,ll)       = (1-zp_sp(subc,ll))*( (1+zp_taus(subc,ll))
                                *zp_gs(subc,ll)
                                + zp_tauc(subc,ll)*calc_tkm('b',subc,ll) )
                                    *(1/calc_tkm('b',subc,ll)) 
                                 + zp_sp(subc,ll)*zp_taup(subc,ll)
                                 + (1+zp_taud(subc,ll))*zp_gd(subc,ll);
    firm_marg_Deltac(subc,ll)        = - calc_v("b",subc,ll)*zp_renth(subc,ll)
                            + calc_vo("b",subc,ll)*firm_marg_Deltamc(subc,ll) 
                            + firm_marg_Deltafc(subc,ll)
                            - calc_pb(subc,ll)*calc_tkm("b",subc,ll)*calc_xB("b",subc,ll);

    firm_marg_costsA(subc,ll)     = zp_renth(subc,ll) * zparam_Hwork(subc)
                                    + zp_dA(subc,ll) ;
    firm_marg_costsB(subc,ll)     = ( zparam_Hwork(subc)  - calc_v('b',subc,ll) )
                                *zp_renth(subc,ll)
                                  + calc_vo('b',subc,ll)*
                                    ( (1-zp_sp(subc,ll))
                                      *( (1+zp_taus(subc,ll))*zp_gs(subc,ll)
                            + zp_tauc(subc,ll)*calc_tkm('b',subc,ll) )
                                      *(1/calc_tkm('b',subc,ll)) 
                                      + zp_sp(subc,ll)*zp_taup(subc,ll)
                                      + (1+zp_taud(subc,ll))*zp_gd(subc,ll)
                                    )
                                  + zp_xbar(subc,ll)*
                                    ( (1+zp_taus(subc,ll))*zp_gs(subc,ll)
                                      + ( (1+zp_taud(subc,ll))*zp_gd(subc,ll)
                                          + zp_tauc(subc,ll)
                                        )*calc_tkm('b',subc,ll)
                                    )
                                  + zp_dB(subc,ll) 
                                  -calc_pb(subc,ll)* calc_tkm('b',subc,ll)
                                     * calc_xB('b',subc,ll);
    firm_chg_marg_costs(subc,ll) = firm_marg_costsA(subc,ll)
                                    - firm_marg_costsB(subc,ll);
    firm_chg_travel(subc,ll)        = calc_v('b',subc,ll) * zp_gs(subc,ll)
                                    - zp_xbar(subc,ll)*zp_gs(subc,ll)
                                      - zp_gs(subc,ll) *calc_tkm('b',subc,ll)*
                                      (calc_xb('a',subc,ll)
                                      -calc_xb('b',subc,ll));

    firm_chg_margprod(subc,ll)  =  zp_beta(subc,ll) * calc_v('b',subc,ll);
    firm_chg_marg_profit(subc,ll)   = firm_chg_margprod(subc,ll) +
                                    firm_chg_marg_costs(subc,ll)
                                    + calc_wb('a',subc,ll)-calc_wb('b',subc,ll);



                                      
* calibrate production function
* benchmark labor demand 
    firm_labor_demand('a',subc,ll) = NLab(subc);
* calibration production parameter
    zp_psi(subc,ll)                 = ( calc_wB('a',subc,ll)+zp_dA(subc,ll)
                                    + Hwork(subc) * zp_renth(subc,ll))
                                    * firm_labor_demand('a',subc,ll);
* includes expected net labor costs considering also change in productivity due to mobile workhours
    firm_avg_expwage('a',subc,ll)     = calc_wb('a',subc,ll)
                                + Hwork(subc) * zp_renth(subc,ll)
                                + zp_dA(subc,ll) ;
                                
    firm_avg_expwage('b',subc,ll)     = calc_wb('a',subc,ll) - ( (zp_adistrib(subc,ll) * power(calc_alpha(subc,ll),2))
                                    /(1-zp_tauw(subc,ll)));
    firm_labor_demand(i,subc,ll)   = max(zp_psi(subc,ll)/firm_avg_expwage(i,subc,ll),0.1);
    
    firm_production(i,subc,ll) =(zp_adistrib(subc,ll) * power(calc_alpha(subc,ll),2))
                                    /(1-zp_tauw(subc,ll));
    firm_production(i,subc,ll)      =
                        zp_psi(subc,ll)*log(firm_labor_demand(i,subc,ll));
    firm_labor_change(subc,ll)     = (firm_labor_demand('b',subc,ll)
                                    - firm_labor_demand('a',subc,ll))
                                /  firm_labor_demand('a',subc,ll);
    firm_margprod(i,subc,ll)       = zp_psi(subc,ll)/firm_labor_demand(i,subc,ll);
    firm_FOC_labordemand(i,subc,ll)= firm_margprod(i,subc,ll) - firm_avg_expwage(i,subc,ll);
    firm_NetProfits(i,subc,ll)     = firm_production(i,subc,ll) - firm_labor_demand(i,subc,ll)
                                    *firm_avg_expwage(i,subc,ll);
    elast_labdemand(subc,ll) = (firm_labor_change(subc,ll) / firm_labor_change(subc,ll))$(firm_labor_change(subc,ll) ne 0);
* effect on office demand per worker

    firm_avg_sqm(i,subc,ll)  = zp_sqmemployee(subc,ll)
                            * ( Hwork(subc) - calc_alpha(subc,ll)*calc_v(i,subc,ll))/Hwork(subc);
    firm_sqm(i,subc,ll) = firm_avg_sqm(i,subc,ll) * firm_labor_demand(i,subc,ll);
    firm_chg_avg_sqm(i,subc,ll) = (firm_avg_sqm('b',subc,ll)
                            -firm_avg_sqm("a",subc,ll) )
                            /firm_avg_sqm("a",subc,ll); 
    firm_chg_sqm(subc,ll) = (firm_sqm("b",subc,ll) -firm_sqm("a",subc,ll) )
                            /firm_sqm("a",subc,ll);

    calc_elastxb(subc,ll) = - calc_pb(subc,ll)* ( 1/
                                (calc_tkm("b",subc,ll)*(calc_uellell(subc,ll)
                                +calc_utxtx(subc,ll))))/
                                (max(calc_xb("b",subc,ll),0.0001));
                                
    calc_avgWTP(subc,ll) = (calc_wB("a",subc,ll) - calc_wB("b",subc,ll))/calc_wB("a",subc,ll)$(calc_wB("a",subc,ll)) + 0$(calc_wB("a",subc,ll)=0);

* ==========================================================================================
* ==========================================================================================
* Print Results to gdx
* ------------------------------------------------------------------------------------------

* set values "0" to eps so that is is printed in gdx file.
    calc_indV(i,subc,ll)$(calc_indV(i,subc,ll)=0) = eps;
    calc_indVx(subc,ll)$(calc_indVx(subc,ll)=0) = eps;
    calc_uvo(subc,ll)$(calc_uvo(subc,ll)=0) = eps;
    calc_uvovo(subc,ll)$(calc_uvovo(subc,ll)=0) = eps; 
    calc_uell(subc,ll)$(calc_uell(subc,ll)=0) = eps;
    calc_uellell(subc,ll)$(calc_uellell(subc,ll)=0) = eps;
    calc_utx(subc,ll)$(calc_utx(subc,ll)=0) = eps;
    calc_utxtx(subc,ll)$(calc_utxtx(subc,ll)=0) = eps; 
    calc_z(i,subc,ll)$(calc_z(i,subc,ll)=0) = eps;
    calc_ell(i,subc,ll)$(calc_ell(i,subc,ll)=0) = eps; 
    calc_vo(i,subc,ll)$(calc_vo(i,subc,ll)=0) = eps;
    calc_effort(i,subc,ll)$(calc_effort(i,subc,ll)=0) = eps;
    calc_v(i,subc,ll)$(calc_v(i,subc,ll)=0) = eps;
    calc_xv(subc,ll)$(calc_xv(subc,ll)=0) = eps;
    calc_hvc(subc,ll)$(calc_hvc(subc,ll)=0) = eps;
    calc_xb(i,subc,ll)$(calc_xb(i,subc,ll)=0) = eps;
    calc_alpha(subc,ll)$(calc_alpha(subc,ll)=0) = eps;
    calc_my(i,subc,ll)$(calc_my(i,subc,ll)=0) = eps;
    calc_my1(subc,ll)$(calc_my1(subc,ll)=0) = eps;
    calc_myc(subc,ll)$(calc_myc(subc,ll)=0) = eps;
    calc_myx(subc,ll)$(calc_myx(subc,ll)=0) = eps;
    calc_myv(subc,ll)$(calc_myv(subc,ll)=0) = eps;
*    calc_my3(subc,ll)$(calc_my3(subc,ll)=0) = eps;
*   calc_my4(subc,ll)$(calc_my4(subc,ll)=0) = eps;
    calc_epsilon(subc,ll)$(calc_epsilon(subc,ll)=0) = eps;
    calc_wB(i,subc,ll)$(calc_wB(i,subc,ll)=0) = eps;
    calc_costA(subc,ll)$(calc_costA(subc,ll)=0) = eps;
    calc_costB(subc,ll)$(calc_costB(subc,ll)=0) = eps;
    calc_wBstar(subc,ll)$(calc_wBstar(subc,ll)=0) = eps;
    calc_tkm(i,subc,ll)$(calc_tkm(i,subc,ll)=0) = eps;
    calc_income(i,subc,ll)$(calc_income(i,subc,ll)=0) = eps;
    calc_avgWTP(subc,ll)$(calc_avgWTP(subc,ll)=0) = eps;
    calc_elastxb(subc,ll)$(calc_elastxb(subc,ll)=0) = eps;
    zindV_z(i,subc,ll)$(zindV_z(i,subc,ll)=0) = eps;
    zindV_ell(i,subc,ll)$(zindV_ell(i,subc,ll)=0) = eps;
    zindV_v(i,subc,ll)$(zindV_v(i,subc,ll)=0) = eps;
    zindV_txB(i,subc,ll)$(zindV_txB(i,subc,ll)=0) = eps;
    calc_pb(subc,ll)$(calc_pb(subc,ll)=0) = eps;
    zp_rho(subc,ll)$(zp_rho(subc,ll)=0) = eps;
    calc_elastxb(subc,ll)$(calc_elastxb(subc,ll)=0) = eps;
    
    
    zemployee_ch_budgetm(i,subc,ll)$(zemployee_ch_budgetm(i,subc,ll)=0) = eps;
    zemployee_ch_budgett(i,subc,ll)$(zemployee_ch_budgett(i,subc,ll)=0) = eps;
    zemployee_tcommute(i,subc,ll)$(zemployee_tcommute(i,subc,ll)=0) = eps;
    zemployee_VOT(i,subc,ll)$(zemployee_VOT(i,subc,ll)=0) = eps;
    zemployee_ut_exputility(subc,ll)$(zemployee_ut_exputility(subc,ll)=0) = eps;
    zemployee_ut_MUt(i,subc,ll)$(zemployee_ut_MUt(i,subc,ll)=0) = eps;
    zemployee_ut_MUxb(i,subc,ll)$(zemployee_ut_MUxb(i,subc,ll)=0) = eps;
    zemployee_ut_MUv(i,subc,ll)$(zemployee_ut_MUv(i,subc,ll)=0) = eps;
    zemployee_ut_MUell(i,subc,ll)$(zemployee_ut_MUell(i,subc,ll)=0) = eps;
    zemployee_ut_MUz(i,subc,ll)$(zemployee_ut_MUz(i,subc,ll)=0) = eps;
    zemployee_ut_MUtxb(i,subc,ll)$(zemployee_ut_MUtxb(i,subc,ll)=0) = eps;
    zemployee_ggen(i,subc,ll)$(zemployee_ggen(i,subc,ll)=0) = eps;
    zemployee_travelcost(i,subc,ll)$(zemployee_travelcost(i,subc,ll) =0) = eps;
    zemployee_traveltime(i,subc,ll)$(zemployee_traveltime(i,subc,ll) =0) = eps;
 
    
    zFOC_xb(i,subc,ll)$(zFOC_xb(i,subc,ll)=0) = eps;
    zFOC_vo(i,subc,ll)$(zFOC_vo(i,subc,ll)=0) = eps;
    zFOC_my1(i,subc,ll)$(zFOC_my1(i,subc,ll)=0) = eps;
    zVOT(i,subc,ll)$(zVOT(i,subc,ll)=0) = eps;
    zVTT(i,subc,ll)$(zVTT(i,subc,ll)=0) = eps; 
    
    zp_xbar(subc,ll)$(zp_xbar(subc,ll)=0) = eps;
    zp_sp(subc,ll)$(zp_sp(subc,ll)=0) = eps;
    zp_dA(subc,ll)$(zp_dA(subc,ll)=0) = eps;
    zp_dB(subc,ll)$(zp_dB(subc,ll)=0) = eps;
    zp_aeff(subc,ll)$(zp_aeff(subc,ll)=0) = eps;
    zp_tauw(subc,ll)$(zp_tauw(subc,ll)=0) = eps;
    zp_tauf(subc,ll)$(zp_tauf(subc,ll)=0) = eps;
    zp_taus(subc,ll)$(zp_taus(subc,ll)=0) = eps;
    zp_taud(subc,ll)$(zp_taud(subc,ll)=0) = eps;
    zp_tauc(subc,ll)$(zp_tauc(subc,ll)=0) = eps;
    zp_taup(subc,ll)$(zp_taup(subc,ll)=0) = eps;
    zp_tauq(subc,ll)$(zp_tauq(subc,ll)=0) = eps;
    zp_rho(subc,ll)$(zp_rho(subc,ll)=0) = eps;
    zp_beta(subc,ll)$(zp_beta(subc,ll)=0) = eps;
    zp_tb(subc,ll)$(zp_tb(subc,ll)=0) = eps;
    zp_y(subc,ll)$(zp_y(subc,ll)=0) = eps;
    zp_renth(subc,ll)$(zp_renth(subc,ll)=0) = eps;
    
    
    firm_NetProfits(i,subc,ll)$(firm_NetProfits(i,subc,ll)=0) = eps;
    firm_chg_travel(subc,ll)$( firm_chg_travel(subc,ll)=0) = eps;
    firm_marg_costsA(subc,ll)$(firm_marg_costsA(subc,ll)=0) = eps;
    firm_marg_costsB(subc,ll)$(firm_marg_costsB(subc,ll)=0) = eps;
    firm_chg_marg_costs(subc,ll)$(firm_chg_marg_costs(subc,ll)=0) = eps;
    firm_marg_Deltafc(subc,ll)$(firm_marg_Deltafc(subc,ll)=0) = eps;
    firm_marg_Deltamc(subc,ll)$(firm_marg_Deltamc(subc,ll)=0) = eps;
    firm_marg_Deltac(subc,ll)$(firm_marg_Deltac(subc,ll)=0) = eps;
    firm_chg_marg_profit(subc,ll)$(firm_chg_marg_profit(subc,ll)=0) = eps;
    firm_chg_margprod(subc,ll)$(firm_chg_margprod(subc,ll)=0) = eps;
    firm_chg_avg_sqm(i,subc,ll)$(firm_chg_avg_sqm(i,subc,ll)=0) = eps;
    firm_FOC_labordemand(i,subc,ll)$(firm_FOC_labordemand(i,subc,ll)=0) = eps;
    elast_labdemand(subc,ll)$(elast_labdemand(subc,ll)=0) = eps;
    
* print to files
$ontext
set
gdxFiles Outputnamen der gdx Files      /MobileWork_Benchmark,MobileWork_tauw,MobileWork_tauf,
                                        MobileWork_tauc,MobileWork_taus,MobileWork_taud,MobileWork_taup,
                                        MobileWork_tauq, MobileWork_rho, MobileWork_gm, MobileWork_gx,
                                        MobileWork_gh, MobileWork_xbar, MobileWork_sp, MobileWork_dA,
                                        MobileWork_dB, MobileWork_rmonth, MobileWork_sqmemployee,
                                        MobileWork_beta, MobileWork_b, MobileWork_y, MobileWork_tb,
                                        MobileWork_tideal, MobileWork_eff, MobileWork_wage,
                                        MobileWork_MonteCarlo
                                        /
s        willkuerliches test set          /s1*s26/
;

parameter iford(gdxFiles)  Dummy fuer ord(iparam) /MobileWork_Benchmark=25,MobileWork_tauw=1,MobileWork_tauf=2,
            MobileWork_tauc=3,MobileWork_taus=4,MobileWork_taud=5,MobileWork_taup=6,
            MobileWork_tauq=7, MobileWork_rho=8, MobileWork_gm=9, MobileWork_gx=10, MobileWork_gh=11,
            MobileWork_xbar=12, MobileWork_sp=13, MobileWork_dA=14, MobileWork_dB=15, MobileWork_rmonth=16,
            MobileWork_sqmemployee=17, MobileWork_beta=18, MobileWork_b=19, MobileWork_y=20, MobileWork_tb=21,
            MobileWork_tideal=22, MobileWork_eff=23, MobileWork_wage=24, MobileWork_MonteCarlo=26/;
            
scalar aa/10/, bb/20/;
$offtext
loop( (s,gdxFiles),
if (ord(s) = ord(iparam),
if(ord(s)=iford(gdxFiles),
   put_utility 'gdxOut' / 'Out' gdxFiles.tl:0;
   execute_unload run_benchmark, 
                run_mc, run_parameters,
                calc_wb, calc_avgWTP, calc_pb, zemployee_ch_budgetm,zemployee_ch_budgett, zemployee_VOT,
                calc_uvo, calc_uvovo, calc_uell, calc_uellell, calc_utx, calc_utxtx,
                i, ctr, ll, calc_z, calc_xb, calc_alpha, calc_income, calc_costA, calc_costB,
                calc_indV, calc_ell, calc_wb, calc_epsilon, calc_wbstar,  calc_indVx,
                calc_v, calc_hvc, calc_vo, calc_myc, calc_myx, calc_myv,
                calc_tkm, calc_elastxb,
                calc_my1, calc_xv, calc_effort, calc_my, calc_Ubar, calc_z,
                calc_elastxb,
*                calc_my3,
*                calc_my4,
                zemployee_ut_exputility, zemployee_ut_MUtxb, zemployee_ut_MUt, 
                zemployee_ut_MUxb, zemployee_ut_MUv, zemployee_ut_MUell, zemployee_ut_MUz,
                zemployee_travelcost,  zemployee_traveltime,zemployee_ggen,
                zemployee_ut_exputility,
                zparam_E_L, zparam_delta0, zparam_delta1, zparam_delta2,
                zparam_delta2, zparam_delta2, zparam_dd2,zparam_delta3,zparam_dd3,
                zparam_peff,
                zparam_Hwork, zparam_Workdays, zp_aeff, zp_tauw, zp_tauf, zp_tauc,
                zp_adistrib, zp_taup, zp_tauq, zp_rho, zp_xbar, zp_timeideal, zp_txideal,
                zp_dA, zp_dB, zp_phi2, zp_varphi,
                zp_rmonth, zp_sqmemployee, zp_beta, zp_gkm, zp_gh, zp_renth,
                zp_sp, zp_avgtime, zp_psi, zp_y, zp_tb,  
                zparam_Hworkmonth, 
                zp_avgproduct, zp_gm, zp_gs, zp_taus, zp_taud, zp_gx, zp_gh, zp_wage,
                zindV_z, zindV_ell, zindV_txb, zindV_v , ZVOT, ZVTT,
                firm_production, firm_margprod, firm_avg_expwage,
                firm_NetProfits, firm_labor_demand,
                firm_labor_change, firm_FOC_labordemand, firm_chg_margprod, firm_chg_travel,
                firm_marg_Deltafc, firm_marg_Deltamc, firm_marg_Deltac, firm_marg_costsA,
                firm_marg_costsB, firm_chg_marg_costs, firm_chg_marg_profit,
                firm_avg_sqm, firm_sqm,
                firm_chg_avg_sqm, firm_chg_sqm, elast_labdemand, tkm0;
);
);
);



);
* ##########################################################################################
*       END SIMULATION and OUTPUT Loop (loop iparam)
* ##########################################################################################
* ##########################################################################################
* ##########################################################################################
