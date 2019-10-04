## Library of functions
import numpy as np
import pdb

def w_mean(data,weigth=None):
    if weigth.any():
        mean=np.sum(data*weigth)/np.sum(weigth)
        dm=np.sqrt(1/np.sum(weigth))
    else:
        mean=np.mean(data)
        dm=np.sqrt(1/len(data))*np.std(data)
    return mean,dm

def chi2(data,dy,model):
    return np.sum(((model - data)**2)/(dy**2))

def linear_regression_AB(x,y,w):
    w=w*np.ones(len(x))
    #pdb.set_trace()
    dw=np.sum(w) * \
            np.sum(w*(x**2)) - \
            (np.sum(w*x))**2
    A=( np.sum(w*(x**2))* \
        np.sum(w*y) - \
        np.sum(w*x) * \
        np.sum(w*x*y) ) / \
        dw
    B=( np.sum(w) * \
        np.sum(w*x*y) - \
        np.sum(w*y) * \
        np.sum(w*x) ) / \
        dw
    dA= np.sqrt(np.sum(((x)**2)*w) / dw)
    dB= np.sqrt(np.sum(w)/dw)
    return A,B,dA,dB

def random_momentum(Ein):
    from real_main import parte1
    mass = parte1.mass
    # momenta are supposed to be gaussianly distributed.
    energy_resolution = 1e-5
    divergence = 1e-5  # all'inizio p è tutto uguale poi il fascio diverge a ventaglio, all'inizio a 0, poi a 1.e-5 o step
    # initial energy
    pout = []
    # prima componente tensore energia impulso è E/c--> c=1
    pout.append(mass + np.random.normal(Ein, energy_resolution*Ein))# energia totale =somma quadratica componenti di p--> c=1
    # momenta
    px_over_p = np.random.normal(0,divergence/np.sqrt(2))
    py_over_p = np.random.normal(0,divergence/np.sqrt(2))
    ps_over_p = np.sqrt(1-px_over_p*px_over_p-py_over_p*py_over_p)
    momentum = np.sqrt(pout[0]*pout[0]-mass*mass) # E^2-p^2c^2=m^2c^4
    pout.append(px_over_p*momentum) #/c con c=1 il momentum
    pout.append(py_over_p*momentum)
    pout.append(ps_over_p*momentum)
    return pout

def random_position():
    # RANDOM_POSITION  Randomly create the initial position of a particle
    # beam profile is supposed to be gaussianly distributed.
    sigmax = 100e-4 # meters
    sigmay = 100e-4
    sigmas = 1e-3
    initial_x_coordinate = 0.015
    initial_y_coordinate = -0.025
    initial_s_coordinate = 0
    xout = []
    xout.append(0) # all particles are created at the same time
    xout.append(np.random.normal(initial_x_coordinate ,sigmax))
    xout.append(np.random.normal(initial_y_coordinate ,sigmay))
    xout.append(np.random.normal(initial_s_coordinate ,sigmas))
    return xout

def RF_cavity(xin, pin):
    from real_main import parte1
    q = parte1.q
    c = parte1.c
    # RF cavity features
    L = 1 #length (meters)
    Ex = 0.001 # Electric field (MV/m)
    Ey = 0.001 # Electric field (MV/m)
    Es = 10.0 # Electric field (MV/m)

    # copy momenta
    pout = []
    xout = []
    pout.append(pin[0])
    pout.append(pin[1])
    pout.append(pin[2])
    pout.append(pin[3])
    xout.append(xin[0])
    xout.append(xin[1])
    xout.append(xin[2])
    xout.append(xin[3])
    # compute relativistic mass
    m2 =  pin[0]*pin[0]-(pin[1]*pin[1]+pin[2]*pin[2]+pin[3]*pin[3])
    # print(m2)
    s0 = 0
    ds = 0.01 # length step (meters)

    while s0 < L:
        # initial velocities
        betax = pout[1]/pout[0] #v/c = p/E=p/mc=p/m
        betay = pout[2]/pout[0]
        betas = pout[3]/pout[0]
        # time interval
        dt = ds/betas  # /light_vel not needed (see momentum's calculation below).. sure?? imho betas'c eliminates this one
        # new coordinates
        xout[0] = xout[0] + dt/c
        xout[1] = xout[1] + betax*dt
        xout[2] = xout[2] + betay*dt
        xout[3] = xout[3] + betas*dt
        # new momenta
        pout[1] = pout[1] + q*0.001*Ex*dt  # units are GeV:  MV-->GV
        pout[2] = pout[2] + q*0.001*Ey*dt
        pout[3] = pout[3] + q*0.001*Es*dt
        pout[0] = np.sqrt(m2+pout[1]*pout[1]+pout[2]*pout[2]+pout[3]*pout[3])
        s0 = s0 + ds
    #print('dt=',dt,'...xout:',xout[0])
    xout[3] = xout[3]-L  # only deltas in s are interesting
    return [xout, pout]

def straight_pipe(xin, pin):
    from real_main import parte1
    c = parte1.c
    xout = []
    # straight pipe length
    L = 20 # meters
    # velocities
    betax = pin[1]/pin[0]
    betay = pin[2]/pin[0]
    betas = pin[3]/pin[0] #s \approx z
    # coordinates
    xout.append(xin[0] + L/betas/c) #dt begin_time +
    dt = xout[0]-xin[0]
    xout.append(xin[1]+betax*c*dt)
    xout.append(xin[2]+betay*c*dt)
    xout.append(betas*c*dt) # only deltas in s are interesting
    # no variation of momenta is expected
    pout = pin # moto rett uniforme
    return [xout, pout]

def compute_dipole(pin): #in teoria non tocco
    from real_main import parte1
    # alpha = L/R = qBL/p
    By = parte1.alpha*pin[3]/parte1.q/parte1.L_dipole
    return By

def dipole(xin, pin, By_in):
    from real_main import parte1
    unit_charge = 1.6e-19;
    # copy momenta
    pout = []
    xout = []
    pout.append(pin[0])
    pout.append(pin[1])
    pout.append(pin[2])
    pout.append(pin[3])
    xout.append(xin[0])
    xout.append(xin[1])
    xout.append(xin[2])
    xout.append(xin[3])
    m2 =  pin[0]*pin[0]-(pin[1]*pin[1]+pin[2]*pin[2]+pin[3]*pin[3])

    s0 = 0
    ds = parte1.dipole_integration_step
    while s0 < parte1.L_dipole:
        #compute the field

        deltax = np.random.normal(0,parte1.disuniformity_dipole/np.sqrt(2))
        deltaz = np.random.normal(0,parte1.disuniformity_dipole/np.sqrt(2))
        compy = np.sqrt(1-deltax*deltax-deltaz*deltaz)
        B = []
        B.append(By_in*deltax) # Megnetic field (T)
        B.append(By_in*compy) # Magnetic field (T)
        B.append(By_in*deltaz) # Magnetic field (T)
        p_vec = []
        p_vec = pout[1:4]

        # initial velocities
        betax = pout[1]/pout[0] #v/c = p/E=p/mc=p/m
        betay = pout[2]/pout[0]
        betas = pout[3]/pout[0]
        beta = [betax,betay,betas]
        #time interval
        dt= ds/betas #  /light_vel   not needed (see momentum' calculation below)
        # new coordinates
        r_vec = []
        r_vec = xout[1:4]
        r_vec = r_vec+dt*np.array(beta)

        # new momenta ###################################
        factor = parte1.q*dt*1.e-9 # units are GeV:  eV-->GeV
        dp_vec = factor*np.cross(beta,B)
        theta = np.linalg.norm(dp_vec)/np.linalg.norm(p_vec)
        p_vec = p_vec + dp_vec
        pout[0] = np.sqrt(m2+np.linalg.norm(p_vec)*np.linalg.norm(p_vec))

        # reprojection on s-direction
        new_s = [-np.sin(theta),0,np.cos(theta)] # local frame with updated s direction
        new_x = np.cross([0,1,0],new_s)/np.linalg.norm(np.cross([0,1,0],new_s))

        pout[3] = np.dot(p_vec,new_s)
        pout[2] = p_vec[1]
        pout[1] = np.dot(p_vec,new_x)

        s_vec = ds*np.array(new_s)
        r_vec = r_vec - s_vec
        xout[0] = xout[0] + dt/parte1.c
        xout[3] = np.dot(r_vec,new_s)
        xout[2] = r_vec[1]
        xout[1] = np.dot(r_vec,new_x)

        s0 = s0 + ds
    return [xout, pout]
