# -*- coding: utf-8 -*-
import assimulo.problem as ap
import numpy as np
import math

class Seven_bar_mechanism(ap.Implicit_Problem):
    """
    A class which describes the squezzer according to
    Hairer, Vol. II, p. 533 ff, see also formula (7.11)
    """
    problem_name = 'Woodpecker w/o friction'
    def __init__(self, index = 3):
        super().__init__()
        self.y0, self.yd0 = self.init_squeezer()
        self.name = self.problem_name
        self.index = index
        self.res = self.f
        
        
    def init_squeezer(self):
        y_1 = np.array([-0.0617138900142764496358948458001,  #  beta
                0.,                                 #  theta
                0.455279819163070380255912382449,   # gamma
                0.222668390165885884674473185609,   # phi
                0.487364979543842550225598953530,   # delta
                -0.222668390165885884674473185609,  # Omega
                1.23054744454982119249735015568])   #epsilon
        lamb = np.array([
            98.5668703962410896057654982170,        # lambda[0]
            -6.12268834425566265503114393122])       # lambda[1]            
        y = np.concatenate((y_1, np.zeros((7,)), lamb, np.zeros((4,))))
        yp = np.concatenate((np.zeros(7,), np.array([
            14222.4439199541138705911625887,        #  betadotdot
            -10666.8329399655854029433719415,       #  Thetadotdot
            0.,0.,0.,0.,0.]), np.zeros((6,))))
        return y, yp
    def f(self,t, y, yp):
        """
        Residual function of the 7-bar mechanism in
        Hairer, Vol. II, p. 533 ff, see also formula (7.11)
        written in residual form
        y,yp vector of dim 20, t scalar
        """
        # Inertia data
        m1,m2,m3,m4,m5,m6,m7=.04325,.00365,.02373,.00706,.07050,.00706,.05498
        i1,i2,i3,i4,i5,i6,i7=2.194e-6,4.410e-7,5.255e-6,5.667e-7,1.169e-5,5.667e-7,1.912e-5
        # Geometry
        xa,ya=-.06934,-.00227
        xb,yb=-0.03635,.03273
        xc,yc=.014,.072
        d,da,e,ea=28.e-3,115.e-4,2.e-2,1421.e-5
        rr,ra=7.e-3,92.e-5
        ss,sa,sb,sc,sd=35.e-3,1874.e-5,1043.e-5,18.e-3,2.e-2
        ta,tb=2308.e-5,916.e-5
        u,ua,ub=4.e-2,1228.e-5,449.e-5
        zf,zt=2.e-2,4.e-2
        fa=1421.e-5
        # Driving torque
        mom=0.033
        # Spring data
        c0=4530.
        lo=0.07785

        # Initial computations and assignments
        beta,theta,gamma,phi,delta,omega,epsilon = y[0:7]
        bep,thp,gap,php,dep,omp,epp = y[7:14]
        lamb = y[14:20]
        sibe,sith,siga,siph,side,siom,siep = np.sin(y[0:7])
        cobe,coth,coga,coph,code,coom,coep = np.cos(y[0:7])
        sibeth = math.sin(beta+theta);cobeth = math.cos(beta+theta)
        siphde = math.sin(phi+delta);cophde = math.cos(phi+delta)
        siomep = math.sin(omega+epsilon);coomep = math.cos(omega+epsilon)
        # Mass matrix
        m = np.zeros((7,7))
        m[0,0] = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 + i2
        m[1,0] = m[0,1] = m2*(da**2-da*rr*coth) + i2
        m[1,1] = m2*da**2 + i2
        m[2,2] = m3*(sa**2+sb**2) + i3
        m[3,3] = m4*(e-ea)**2 + i4
        m[4,3] = m[3,4] = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
        m[4,4] = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2)+ i4 + i5
        m[5,5] = m6*(zf-fa)**2 + i6
        m[6,5] = m[5,6] = m6*((zf-fa)**2-u*(zf-fa)*siom) + i6
        m[6,6] = m6*((zf-fa)**2-2*u*(zf-fa)*siom+u**2) + m7*(ua**2+ub**2)+ i6 + i7

        # Applied forces

        xd = sd*coga + sc*siga + xb
        yd = sd*siga - sc*coga + yb
        lang  = math.sqrt((xd-xc)**2 + (yd-yc)**2)
        force = - c0 * (lang - lo)/lang
        fx = force * (xd-xc)
        fy = force * (yd-yc)
        ff = np.array([
        mom - m2*da*rr*thp*(thp+2*bep)*sith,    
        m2*da*rr*bep**2*sith,
        fx*(sc*coga - sd*siga) + fy*(sd*coga + sc*siga),
        m4*zt*(e-ea)*dep**2*coph,
        - m4*zt*(e-ea)*php*(php+2*dep)*coph,
        - m6*u*(zf-fa)*epp**2*coom,
        m6*u*(zf-fa)*omp*(omp+2*epp)*coom])

        #  constraint matrix  G

        gp = np.zeros((6,7))
    
        gp[0,0] = - rr*sibe + d*sibeth
        gp[0,1] = d*sibeth
        gp[0,2] = - ss*coga
        gp[1,0] = rr*cobe - d*cobeth
        gp[1,1] = - d*cobeth
        gp[1,2] = - ss*siga
        gp[2,0] = - rr*sibe + d*sibeth
        gp[2,1] = d*sibeth
        gp[2,3] = - e*cophde
        gp[2,4] = - e*cophde + zt*side
        gp[3,0] = rr*cobe - d*cobeth
        gp[3,1] = - d*cobeth
        gp[3,3] = - e*siphde
        gp[3,4] = - e*siphde - zt*code
        gp[4,0] = - rr*sibe + d*sibeth
        gp[4,1] = d*sibeth
        gp[4,5] = zf*siomep
        gp[4,6] = zf*siomep - u*coep
        gp[5,0] = rr*cobe - d*cobeth
        gp[5,1] = - d*cobeth
        gp[5,5] = - zf*coomep
        gp[5,6] = - zf*coomep - u*siep

        #     Index-3 constraint
        g = np.zeros((6,))
        g[0] = rr*cobe - d*cobeth - ss*siga - xb
        g[1] = rr*sibe - d*sibeth + ss*coga - yb
        g[2] = rr*cobe - d*cobeth - e*siphde - zt*code - xa
        g[3] = rr*sibe - d*sibeth + e*cophde - zt*side - ya
        g[4] = rr*cobe - d*cobeth - zf*coomep - u*siep - xa
        g[5] = rr*sibe - d*sibeth - zf*siomep + u*coep - ya
        
        # Construction of g_qq
        g_qq = np.zeros(6)
        g_qq[0] = -rr*cobe*y[7]**2+d*cobeth*(y[7]+y[8])**2+ss*siga*y[9]**2
        g_qq[1] = -rr*sibe*y[7]**2+d*sibeth*(y[7]+y[8])**2-ss*coga*y[9]**2
        g_qq[2] = -rr*cobe*y[7]**2+d*cobeth*(y[7]+y[8])**2+e*siphde*(y[10]+y[11])**2+zt*code*y[11]**2
        g_qq[3] = -rr*sibe*y[7]**2+d*sibeth*(y[7]+y[8])**2-e*cophde*(y[10]+y[11])**2+zt*side*y[11]**2
        g_qq[4] = -rr*cobe*y[7]**2+d*cobeth*(y[7]+y[8])**2+zf*coomep*(y[12]+y[13])**2+u*siep*y[13]**2
        g_qq[5] = -rr*sibe*y[7]**2+d*sibeth*(y[7]+y[8])**2+zf*siomep*(y[12]+y[13])**2-u*coep*y[13]**2

        #     Construction of the residual
    
        res_1 = yp[0:7] - y[7:14]
        res_2 = np.dot(m, yp[7:14])- ff[0:7] + np.dot(gp.T,lamb)
        if self.index == 3:
            res_3 = g
        elif self.index == 2:
            res_3 = np.dot(gp, y[7:14])
        else:
            raise Exception('Only index 2 or 3 index formulations available')
        
        return np.concatenate((res_1, res_2, res_3))
    

    
    
    


    
