import numpy as np
from uncertainties import ufloat
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from functions import *
global mass, n, q, c

def parte1():
    ### initialize_particle
    parte1.q = 1 #1.6e-19
    parte1.c = 2.998e8
    parte1.mass = 0.938 #rest mass of the particle GeV/c^2
    parte1.n = 5 #number of particles
    n = parte1.n

    ### initialize_dipoles
    parte1.N_dipoles = 2
    parte1.L_dipole = 2
    parte1.dipole_integration_step = 0.05 #meters, inizio con tanto, se dipolo lungo 1 metro, lo dividio in 100 passi ca
    parte1.alpha = 2*np.pi/parte1.N_dipoles #every dipole should bend the beam of this angle
    parte1.disuniformity_dipole = 1e-6

    ### initialize_quadrupoles
    L_quadrupole = .5 # length (meters)
    quadrupole_integration_step = 0.005 #meters
    disuniformity_quadrupole = 0 #1.e-6 #not used by now

    # starting_point
    X = []
    P = []
    fig1 = plt.figure()
    for i in range(0, n):
        pstart = random_momentum(450) #450 gev iniziali
        rstart = random_position()
        X.append(rstart)
        P.append(pstart)
        print('P iniziale della',i+1,'particella:')
        print(P[i])

    for idx in range(0,1000):
        #print(idx)
        for i in range(0,n):
            [X[i], P[i]] = RF_cavity(X[i],P[i])
            [X[i], P[i]] = straight_pipe(X[i],P[i])
            By = compute_dipole(P[i])
            dipole(X[i], P[i], By)
            [X[i], P[i]] = straight_pipe(X[i],P[i])
            colours = 'byrcmgk'
            By = compute_dipole(P[i])
            print('By=',By)
            dipole(X[i], P[i], By)

            plt.plot(X[i][1],P[i][2], colours[i]+'o')
            plt.grid()
            axes = plt.gca()
            axes.set_xlabel(idx)
            axes.set_xlim([-0.1,0.1])
            axes.set_ylim([-0.1,0.1])

        fig1.show()
        plt.pause(0.001)
        fig1.clear()
    print(X[1],P[1])


# for idx=1:1000
#     for i=1:5
#     [X{i},P{i}] = RF_cavity(X{i},P{i});
#     [X{i},P{i}] = straight_pipe(X{i},P{i});
#     pp = plot(X{i}(2),X{i}(3),strcat(colours(i),'o'));
#
#     g_in = compute_quadrupole(P{i});
#     [X{i},P{i}] = quadrupole(X{i},P{i},1 ,g_in);
#
#
#     [X{i},P{i}] = straight_pipe(X{i},P{i});
#     [X{i},P{i}] = quadrupole(X{i},P{i},0 ,g_in);
#
#     By = compute_dipole(P{i});
#     [X{i},P{i}] = dipole(X{i},P{i},By);
#
#     [X{i},P{i}] = straight_pipe(X{i},P{i});
#
#     [X{i},P{i}] = dipole(X{i},P{i},By);
#
#     grid on;
#     xlim([-0.1,0.1]);
#     ylim([-0.1,0.1]);
#     hold on
#     end
#     pause(0.001);
#     %hold off; %commentare se si vuole vedere la traiettoria
# end
# for i=1:5
#     fprintf(['P finale della ' int2str(i) ' particella:']);
#     disp(P{i});
# end
# hold off


    '''
    def import_misure(file_name):
        with open(file_name) as csvfile:
            data=[]
            reader=csv.reader(csvfile,quoting=csv.QUOTE_NONNUMERIC)
            for i,row in enumerate(reader):
                if row and (not isinstance(row[1],str)):
                    data.append(row)
            data=np.array(data)
            return data[:,0],data[:,1],data[:,2]
    def sat_vapo_pres(t):
        return (10**(a-(b/(c+t))))*133.322

    def covar(t, w):
        return (np.sum(t*w)/((np.sum(w)*np.sum(w*t*t))-(np.sum(w*t))**2))

    ## UNCERT
    dh_m=1e-3/np.sqrt(6)
    dh_a=1e-3/np.sqrt(6)
    dtemp=1e-2/np.sqrt(12)

    ## DATA
    time_m,h_m,temp=import_misure('data/misure_lun_14_5_18.csv')
    h_m = 95-h_m
    h_m = 1e-2*h_m
    # temp+=ZERO_C
    time_a,h_a,pres_a=import_misure('data/pr_atm_lun_14_5_18.csv')
    #TODO: set 72.2 for height of mercury
    h_a=72.4*np.ones(h_a.shape)
    h_a*=1e-2
    pres_a*=1e2
    # ##calculating atmospheric pression at measured time
    p_a=RHO_HG*G_ACC*h_a+pres_a
    dp_a=RHO_HG*G_ACC*dh_a
    p_a_mt=np.empty(h_m.shape)
    for i in range(len(p_a_mt)):###
        #i_b is the index in time_a for the time before time_m[i]
        i_b=0
        while(not (time_m[i] < time_a[i_b])):
            i_b+=1
        i_b-=1
        p_a_mt[i]=p_a[i_b]+(p_a[i_b+1]-p_a[i_b])/(time_a[i_b+1]-time_a[i_b])*(time_m[i]-time_a[i_b])
    ##nice perameters
    i_zero=np.argmin(temp)
    dh_m_trasf = 3.60*dtemp
    h_m_wei=1/(np.sqrt((dh_m**2)+(dh_m_trasf**2)))

    p_m=p_a_mt-RHO_H2O*G_ACC*h_m
    dp_m=np.sqrt(dp_a**2+(RHO_H2O*G_ACC*dh_m)**2)
    dp_m_trasf = 354*dtemp
    p_m_wei=1/((dp_m**2)+(dp_m_trasf**2))
    p_0=p_m[i_zero]
    alpha=1/ZERO_C


    #########################################################################
    #             ANALYSIS
    #########################################################################
    i_min=np.argmin(temp)
    A1,B1,dA1,dB1=linear_regression_AB(temp[:],h_m[:],h_m_wei)
    print('Linear regression dz e t: h =',ufloat(A1,dA1),'+',ufloat(B1,dB1),'* T')
    chi2_ht1=chi2(h_m,dh_m,A1+B1*temp)
    print('Chi2 dz e T: ',chi2_ht1)

    A2,B2,dA2,dB2=linear_regression_AB(temp[:],p_m[:],p_m_wei)
    print('Linear regression P e t: P =',ufloat(A2,dA2),'+',ufloat(B2,dB2),'* T')
    chi2_ht2=chi2(p_m,dp_m,A2+B2*temp)
    print('Chi2 dz e T: ',chi2_ht2)

    # print('Considero solo i valori con T>D_TEMP')
    temp2=[]
    h_m2=[]
    p_m2=[]
    for i,t in enumerate(temp):
        if t > D_TEMP:
            temp2.append(t)
            h_m2.append(h_m[i])
            p_m2.append(p_m[i])
    temp2=np.array(temp2)
    h_m2=np.array(h_m2)
    p_m2=np.array(p_m2)
    i_min2=np.argmin(temp2)

    temp3=[]
    h_m3=[]
    p_m3=[]
    for i,t in enumerate(temp):
        if t < D_TEMP+0.5:
            temp3.append(t)
            h_m3.append(h_m[i])
            p_m3.append(p_m[i])
    temp3=np.array(temp3)
    h_m3=np.array(h_m3)
    p_m3=np.array(p_m3)
    i_min3=np.argmin(temp3)


    A3,B3,dA3,dB3=linear_regression_AB(temp2[:],p_m2[:],p_m_wei)
    print('Linear regression sopra: P =',A3,'+/-',dA3,'+',B3, '+/-', dB3,'* T')
    chi2_ht3=chi2(p_m2[:],dp_m,A3+B3*temp2[:])
    print('Chi 3: ',chi2_ht3)

    A4,B4,dA4,dB4=linear_regression_AB(temp3[:],p_m3[:],p_m_wei)
    print('Linear regression sotto: h =',ufloat(A4,dA4),'+',ufloat(B4,dB4),'* T')
    chi2_ht4=chi2(p_m3[:],dp_m,A4+B4*temp3[:])
    print('Chi 4: ',chi2_ht4)


    # A3,B3,dA3,dB3=linear_regression_AB(temp2,p_m2,p_m_wei)
    # print('Linear regression total: P =',ufloat(A3,dA3),'+',ufloat(B3,dB3),'* T')
    # chi2_pt3=chi2(p_m2,dp_m,A3+B3*temp2)
    # print('Chi2 3: ',chi2_pt3)
    #
    # A31,B31,dA31,dB31=linear_regression_AB(temp2[:i_min2],p_m2[:i_min2],p_m_wei)
    # print('Linear regression scendendo: P =',ufloat(A31,dA31),'+',ufloat(B31,dB31),'* T')
    # chi2_pt31=chi2(p_m2[:i_min2],dp_m,A31+B31*temp2[:i_min2])
    # print('Chi2 31: ',chi2_pt31)
    #
    # A32,B32,dA32,dB32=linear_regression_AB(temp2[i_min2:],p_m2[i_min2:],p_m_wei)
    # print('Linear regression salendo: P =',ufloat(A32,dA32),'+',ufloat(B32,dB32),'* T')
    # chi2_pt32=chi2(p_m2[i_min2:],dp_m,A32+B32*temp2[i_min2:])
    # print('Chi2 32: ',chi2_pt32)
    #
    # #TODO:redo when chi2 works
    #
    # #let's try to use under dew point data
    a=8.07131
    b=1730.63
    c=233.426
    # # #abs_hum is in kg*m^-3
    # # abs_hum=2.16679e-3*sat_vapo_pres(20+ZERO_C)/(20+ZERO_C)*RH
    # # #percentage of water in air by mass
    # # rho_air_20C=1.2041
    # # water_in_air=abs_hum/rho_air_20C
    # #correcting pressure adding the help of water that condensed
    t_comoda = 20
    p_m_cor=p_m.copy()
    for i in range(len(p_m_cor)):
        if temp[i] < 12+0.5:
            # p_m_cor[i]+=water_in_air*(sat_vapo_pres(D_TEMP)+(temp[i]-D_TEMP)*(sat_vapo_pres(D_TEMP)/temp[i])-sat_vapo_pres(temp[i]))
            #TODO:it works like that but it doesn't make sense at all
            # p_m_cor[i]+=RH*(sat_vapo_pres(D_TEMP)+(temp[i]-D_TEMP)*(sat_vapo_pres(D_TEMP)/temp[i])-sat_vapo_pres(temp[i]))
            # sottraggo al valore della retta in temp[i] la funzione esponenziale in temp[i]
            p_m_cor[i]+=((.66*sat_vapo_pres(t_comoda)/(t_comoda*273.15))*temp[i]+273.15*(.66*sat_vapo_pres(t_comoda)/(t_comoda+273.15)))-sat_vapo_pres(temp[i])
            # print(temp[i], ((.66*sat_vapo_pres(t_comoda)/(t_comoda*273.15))*temp[i]+273.15*(.66*sat_vapo_pres(t_comoda)/(t_comoda+273.15)))-sat_vapo_pres(temp[i]))

    A5,B5,dA5,dB5=linear_regression_AB(temp,p_m_cor,p_m_wei)
    print('Linear regression corr: P =',ufloat(A5,dA5),'+',ufloat(B5,dB5),'* T')
    chi2_pt5=chi2(p_m_cor,dp_m,A5+B5*temp)
    print('Chi2 corretto pressioni:',chi2_pt5)

    A6,B6,dA6,dB6=linear_regression_AB(temp2[:9],p_m2[:9],p_m_wei)
    A7,B7,dA7,dB7=linear_regression_AB(temp2[10:],p_m2[10:],p_m_wei)

    w_cov = 1/(dp_m)**2
    zero_guess1=-A3/B3
    dzero_guess1=np.sqrt((dA3/B3)**2+(dB3*A3/(B3**2))**2-(2*A3*covar(temp2, w_cov)/B3**3))
    print('zero guess 1:',zero_guess1,'+/-', dzero_guess1)

    w_cov = 1/(dp_m)**2
    zero_guess2=-A5/B5
    dzero_guess2=np.sqrt((dA5/B5)**2+(dB5*A5/(B5**2))**2-(2*A5*covar(temp, w_cov)/B5**3)**2)
    dzero_guess3=dA3/B3+dB3*A3/B3**2
    dzero_guess4=np.sqrt((dA3/B3)**2+(dB3*A3/B3**2)**2)


    print('zero guess 1:',zero_guess2,'+/-', dzero_guess2)
    print('con schwartz:',zero_guess2,'+/-', dzero_guess3)
    print('senza covarianza:',zero_guess2,'+/-', dzero_guess4)


    # zero_guess2=-A4/B4
    # dzero_guess2=np.sqrt((dA4/B4)**2+(A4/B4**2*dB4)**2)
    # print('zero guess 2:',ufloat(zero_guess2,dzero_guess2))
    #
    # print('dp_a:',dp_a)
    ############################################################################
    #             PLOTS
    ############################################################################
    #PLOT1 - GRAFICO DI H IN FUNZIONE DI T
    fig1=plt.figure(figsize=DOUBLE_FIGSIZE)
    ax11=fig1.add_subplot(1,1,1)
    ax11.errorbar(temp[:i_zero],h_m[:i_zero]*1e2,yerr=dh_m,xerr=dtemp,fmt='b.',label='Dati andata')
    ax11.errorbar(temp[i_zero:],h_m[i_zero:]*1e2,yerr=dh_m,xerr=dtemp,fmt='r.',label='Dati ritorno')
    # ax12=ax11.twinx()
    # ax12.set_ylabel('Pgas [Pa]')
    # ax12.errorbar(temp,p_m,yerr=dh_m,xerr=dtemp,fmt='g.',label='Temp salendo')
    ax11.set_xlabel('T [°C]')
    ax11.set_ylabel('$\Delta z\ [cm]$')
    ax11.grid()
    ax11_zoom=fig1.add_axes((.23,.23,.2,.2))
    ax11_zoom.errorbar(temp[11],h_m[11]*1e2,yerr=dh_m,xerr=dtemp, fmt = 'b')
    ax11_zoom.set_xlabel('T [°C]')
    ax11_zoom.set_ylabel('$\Delta z\ [cm]$')
    ax11_zoom.grid()
    mark_inset(ax11, ax11_zoom, loc1=2, loc2=4, fc="none", ec="0.5")
    fig1.suptitle('$\Delta z$ in funzione della temperatura',fontsize=18)
    legend1 = ax11.legend(loc='upper right', shadow=True, prop={'size': 17})
    legend1.get_frame().set_facecolor('#FFFFFF')
    fig1.savefig('fig1.png', transparent=False, dpi=180, )


    # #PLOT2 - GRAFICO RESIDUI DI H IN FUNZIONE DI T SOPRA D-POINT
    fig2=plt.figure(figsize=DOUBLE_FIGSIZE)
    fig2.suptitle('Residuo normalizzato in funzione della temperatura',fontsize=17)
    ax21=fig2.add_subplot(1,1,1)
    ax21.errorbar(temp[:i_zero],(h_m[:i_zero]-(A1+B1*temp[:i_zero]))/dh_m,yerr=dh_m,xerr=dtemp,fmt='b.',label='Dati andata')
    ax21.errorbar(temp[i_zero:],(h_m[i_zero:]-(A1+B1*temp[i_zero:]))/dh_m,yerr=dh_m,xerr=dtemp,fmt='r.',label='Dati ritorno')
    ax21.axhline(y=0)
    ax21.grid()
    ax21.set_xlabel('T [°C]')
    ax21.set_ylabel('$R/\sigma [\Delta z]$')
    legend2 = ax21.legend(loc='lower left', shadow=True, prop={'size': 15})
    legend2.get_frame().set_facecolor('#FFFFFF')
    fig2.savefig('fig2.png', transparent=False, dpi=180, )

    #PLOT3 - PRESSIONI NON CORRETTE SOPRA D-POINT SCENDENDO
    fig3=plt.figure(figsize=DOUBLE_FIGSIZE)
    ax31=fig3.add_subplot(1,1,1)
    ax31.errorbar(temp[:i_zero],p_m_cor[:i_zero],yerr=dp_m,xerr=dtemp,fmt='b.',label='Dati andata')
    ax31.errorbar(temp[i_zero:],p_m_cor[i_zero:],yerr=dp_m,xerr=dtemp,fmt='r.',label='Dati ritorno')
    ax31.set_xlabel('T [°C]')
    ax31.set_ylabel('P [Pa]')
    fig3.suptitle('Pressione in funzione della temperatura',fontsize=18)
    ax31.grid()
    ax31_zoom=fig3.add_axes((.65,.23,.2,.2))
    ax31_zoom.errorbar(temp[11],p_m_cor[11],yerr=dp_m,xerr=dtemp, fmt = 'b')
    ax31_zoom.set_xlabel('T [°C]')
    ax31_zoom.set_ylabel('P [Pa]')
    ax31_zoom.grid()
    mark_inset(ax31, ax31_zoom, loc1=1, loc2=3, fc="none", ec="0.5")
    legend3 = ax31.legend(loc='upper left', shadow=True, prop={'size': 17})
    legend3.get_frame().set_facecolor('#FFFFFF')
    fig3.savefig('fig3.png', transparent=False, dpi=180, )

    #GRAFICO ANTOINE
    fig4=plt.figure(figsize=DOUBLE_FIGSIZE)
    ax41=fig4.add_subplot(1,1,1)
    x4a = np.linspace(-2, 12.35)
    y4a = sat_vapo_pres(x4a)
    x4b = np.linspace(-2, 12.35)
    y4b = ((.66*sat_vapo_pres(t_comoda)/(t_comoda*273.15))*x4b+273.15*(.66*sat_vapo_pres(t_comoda)/(t_comoda+273.15)))
    x4c = np.linspace(12.35, 23)
    y4c = ((.66*sat_vapo_pres(t_comoda)/(t_comoda*273.15))*x4c+273.15*(.66*sat_vapo_pres(t_comoda)/(t_comoda+273.15)))
    x4d = np.linspace(12.35, 23)
    y4d = sat_vapo_pres(x4d)
    ax41.grid()
    ax41.set_xlabel('T [°C]')
    ax41.set_ylabel('P [Pa]')
    fig4.suptitle('Pressione di vapore saturo - H2O',fontsize=20)
    ax41.plot(x4a, y4a, c='r', lw=2.5, label = '$P_{H2O}$ effettiva')
    ax41.plot(x4d, y4d, c='b', ls='--', lw=2, label = 'Pressione di vapor saturo')
    ax41.plot(x4b, y4b, c='g', ls='--', lw=2, label = 'Modello per gas perfetto')
    ax41.plot(x4c, y4c, c='r', lw=2.5)
    legend4 = ax41.legend(loc='upper left', shadow=True, prop={'size': 14})
    legend4.get_frame().set_facecolor('#FFFFFF')
    fig4.savefig('fig4.png', transparent=False, dpi=180, )

    # #GRAFICI RESIDUI
    # fig5=plt.figure(figsize=DOUBLE_FIGSIZE)
    # x1 = np.linspace(0.0, 5.0)
    # x2 = np.linspace(0.0, 2.0)
    # y1 = np.cos(2 * np.pi * x1) * np.exp(-x1)
    # y2 = np.cos(2 * np.pi * x2)
    # plt.subplot(1, 2, 1)
    # plt.plot(x1, y1, 'ko-')
    # plt.title('A tale of 2 subplots')
    # plt.ylabel('Damped oscillation')
    # plt.subplot(1, 2, 2)
    # plt.plot(x2, y2, 'r.-')
    # plt.xlabel('time (s)')
    # plt.ylabel('Undamped')
    # plt.tight_layout()



    # #PLOT4 - PRESSIONI NON CORRETTE SOPRA D-POINT SALENDO
    # fig4=plt.figure(figsize=DOUBLE_FIGSIZE)
    # fig4.suptitle('Dipendenza di P da T salendo',fontsize=16)
    # ax41=fig4.add_subplot(1,1,1)
    # ax41.errorbar(temp2[i_min2:],p_m2[i_min2:]-A32-B32*temp2[i_min2:],yerr=dp_m,xerr=dtemp,fmt='b.',label='Temp scendendo')
    # ax41.axhline(y=0)
    # ax41.set_xlabel('T [K]')
    # ax41.set_ylabel('P [Pa]')
    # #
    # #PLOT5 - PRESSIONI 40 MISURE CON CORREZIONE ANTOINE
    # fig5=plt.figure(figsize=DOUBLE_FIGSIZE)
    # fig5.suptitle('Pressioni corrette tenendo conto della pressione di vapore saturo',fontsize=16)
    # ax51=fig5.add_subplot(1,1,1)
    # ax51.errorbar(temp,p_m_cor,xerr=dtemp,yerr=dp_m,fmt='.',label='Pres corr')
    # ax51.errorbar(temp,p_m,xerr=dtemp,yerr=dp_m,fmt='.',label='Pres orig')
    # ax51.set_xlabel('T [K]')
    # ax51.set_ylabel('P [Pa]')
    #
    # #PLOT6 - RESIDUI PRESSIONI CORRETTE
    # fig6=plt.figure(figsize=DOUBLE_FIGSIZE)
    # fig6.suptitle('Residui pressioni corrette',fontsize=16)
    # ax61=fig6.add_subplot(1,1,1)
    # ax61.errorbar(temp[i_zero:],p_m_cor[i_zero:]-A4-B4*temp[i_zero:],xerr=dtemp,yerr=dp_m,fmt='b.',label='Residui pres corr')
    # ax61.errorbar(temp[:i_zero],p_m_cor[:i_zero]-A4-B4*temp[:i_zero],xerr=dtemp,yerr=dp_m,fmt='r.',label='Residui pres corr')
    # ax61.plot(temp[i_zero:],p_m_cor[i_zero:]-A4-B4*temp[i_zero:],'-')
    # ax61.plot(temp[:i_zero],p_m_cor[:i_zero]-A4-B4*temp[:i_zero],'-')
    # ax61.axhline(y=0,color='r')
    # ax61.axvline(x=D_TEMP,color='g')
    # ax61.set_xlabel('T [K]')
    # ax61.set_ylabel('P [Pa]')
    #
    plt.show() '''
