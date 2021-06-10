import numpy as np
import matplotlib.pyplot as plt


class Filament:

    def __init__(self, z_fil, z_cp, w, V_inf, G_A, y_A, phi_p, phi_G):

        # Store vars
        self.z_fil = z_fil
        self.z_cp = z_cp
        self.w = w
        self.V_inf = V_inf
        self.G_A = G_A
        self.y_A = y_A
        self.phi_p = phi_p
        self.phi_G = phi_G

        self.dz = z_cp-z_fil
        self.w_y = w*y_A/V_inf

    
    def S_p(self, x):
        return np.sin(self.w*x/self.V_inf+phi_p)

    
    def C_p(self, x):
        return np.cos(self.w*x/self.V_inf+phi_p)

    
    def S_G(self, x):
        return np.sin(self.w*x/self.V_inf+phi_G)

    
    def C_G(self, x):
        return np.cos(self.w*x/self.V_inf+phi_G)


    def dV(self, x):
        # Calculates the differential element of induced velocity from a station x on the filament

        # Get sinusoids
        S_G = self.S_G(x)
        C_p = self.C_p(x)
        S_p = self.S_p(x)

        # Initialize result vector
        dV = np.zeros(3)

        # Compute radial distance
        r = (x**2+self.y_A**2*S_p**2+self.dz**2)**1.5

        # Compute cross product
        dV[0] = self.dz*self.w_y*C_p
        dV[1] = -self.dz
        dV[2] = -self.y_A*S_p+x*self.w_y*C_p

        return self.G_A*S_G*dV/(4.0*np.pi*r)


    def dV_alt(self, theta):
        # Calculates the differential element of induced velocity from a station theta on the filament assuming a straight filament

        # Calculate params
        h = self.dz

        if theta == 0.0:
            return 0.0
        else:
            return -self.G_A*np.sin(theta)*np.sin(self.w*h/(self.V_inf*np.tan(theta))+self.phi_G)/(4.0*np.pi*h)


def gauss_quad_sixth_integration_vec(f, a, b):
    # Uses sixth-order Gaussian quadrature to compute the integral of a function

    # Formula constants
    c = np.array([0.1713245, 0.3607616, 0.4679139, 0.4679139, 0.3607616, 0.1713245])
    x = np.array([-0.932469514, -0.661209386, -0.238619186, 0.238619186, 0.661209386, 0.932469514])

    # Transform to integration domain
    x_t = 0.5*x*(b-a)+a+1

    # Get function values
    vals = np.zeros((len(x), 3))
    for i, x_i in enumerate(x):
        vals[i] = f(x_t[i])

    # Sum and rescale
    return 0.5*(b-a)*np.sum(c[:,np.newaxis]*vals, axis=0)


def boole_integration(f, a, b, N):
    # Uses Boole's rule to compute the integral of f from a to b using N points (N must be N = 5*m+1)

    # Get function values
    x = np.linspace(a, b, N)
    vals = np.zeros(N)
    for i in range(N):
        vals[i] = f(x[i])

    # Get summation
    I = 0.0
    for i in range(0, N-1, 5):
        I += (x[i+5]-x[i])*(19.0*vals[i]+75.0*vals[i+1]+50.0*vals[i+2]+50.0*vals[i+3]+75.0*vals[i+4]+19.0*vals[i+5])

    return I/288.0


def boole_integration_vec(f, a, b, N):
    # Uses Boole's rule to compute the integral of f from a to b using N points (N must be N = 5*m+1)

    # Get function values
    x = np.linspace(a, b, N)
    vals = np.zeros((N, 3))
    for i in range(N):
        vals[i] = f(x[i])

    # Get summation
    I = np.zeros(3)
    for i in range(0, N-1, 5):
        I += (x[i+5]-x[i])*(19.0*vals[i]+75.0*vals[i+1]+50.0*vals[i+2]+50.0*vals[i+3]+75.0*vals[i+4]+19.0*vals[i+5])

    return I/288.0


if __name__=="__main__":

    # Options
    G_A = 1.0
    z0 = 1.0
    z = 0.0
    w = 1.0
    V_inf = 10.0
    y_A = 0.0
    phi_p = 0.0
    phi_G = 0.0

    # Initialize filament
    fil = Filament(z, z0, w, V_inf, G_A, y_A, phi_p, phi_G)

    # Compare integrands
    V = boole_integration_vec(fil.dV, 0.0, 1000.0, 100001)
    print("Original: Vy = {0}".format(V[1]))
    V = boole_integration(fil.dV_alt, 0.5*np.pi, 0.0, 100001)
    print("Simplified: Vy = {0}".format(V))

    ## Loop through limits of integration
    #bs = np.logspace(0, 4, 50)
    #Vs = np.zeros((50,3))
    #for i, b in enumerate(bs):
    #    Vs[i] = gauss_quad_sixth_integration_vec(fil.dV, 0, b)

    ## Plot results
    #plt.figure()
    #plt.plot(bs, Vs[:,0], label='Vx')
    #plt.plot(bs, Vs[:,1], label='Vz')
    #plt.plot(bs, Vs[:,2], label='Vy')
    #plt.title("Gaussian Quadrature")
    #plt.legend()
    #plt.xlabel('b')
    #plt.ylabel('Induced Vel')
    #plt.xscale('log')
    #plt.show()

    ## Loop through limits of integration
    #dx_des = 0.001
    #bs = np.logspace(0, 3, 10)
    #Vs = np.zeros((10,3))
    #for i, b in enumerate(bs):
    #    print(b)
    #    N = b/dx_des
    #    N = int(N//5)
    #    N = N*5+1
    #    Vs[i] = boole_integration_vec(fil.dV, 0, b, N)

    ## Plot results
    #plt.figure()
    #plt.plot(bs, Vs[:,0], label='Vx')
    #plt.plot(bs, Vs[:,1], label='Vy')
    #plt.plot(bs, Vs[:,2], label='Vz')
    #plt.title("Boole's Rule")
    #plt.legend()
    #plt.xlabel('b')
    #plt.ylabel('Induced Vel')
    #plt.xscale('log')
    #plt.show()

    #print(Vs[-1,:])
