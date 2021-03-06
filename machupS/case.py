import json

import numpy as np
import matplotlib.pyplot as plt

from machupS.wing import Wing, vec_norm, vec_cross, vec_inner


class Case:
    """Defines a steady-periodic wing simulation."""

    def __init__(self, **kwargs):

        # Get kwargs
        wing_kwargs = kwargs.get("wing")

        state_kwargs = kwargs.get("state")
        self._V_inf = state_kwargs.get("V_inf")
        self._w = state_kwargs.get("frequency")
        self._z_A = state_kwargs.get("plunge_amplitude")
        self._phi_z = np.radians(state_kwargs.get("plunge_phase"))
        self._a_A = np.radians(state_kwargs.get("pitch_amplitude"))
        self._phi_a = np.radians(state_kwargs.get("pitch_phase"))
        self._a_mean = np.radians(state_kwargs.get("mean_alpha"))
        
        case_options = kwargs.get("case_options")
        self._max_iter = case_options.get("max_iter")
        planar = case_options.get("planar")
        self._dt = case_options.get("dt")

        # Determine initial velocity
        V_init = np.zeros(3)
        V_init[0] = self._V_inf
        V_init[2] = 2.0*np.pi*self._z_A*self._w/self._V_inf*np.cos(self._phi_z)

        # Initialize wing
        self.wing = Wing(**wing_kwargs, max_iter=self._max_iter, planar=planar, V_init=V_init, V_inf=self._V_inf)


    def run(self, results_file='machups_results.csv', verbose=False):
        """Runs the case as specified."""

        # Initialize storage
        V_sig = np.zeros((self.wing.N, 3))
        self.gamma = np.zeros((self._max_iter, self.wing.N))

        # Initialize command line output
        if verbose:
            print("{0:<20}{1:<20}{2:<20}{3:<20}{4:<20}".format("Iteration", "Pitch Angle [deg]", "Plunge Vel [ft/s]", "CL", "CDi"))
            print("".join(["-"]*100))

        # Open results file
        with open(results_file, 'w') as results_file_handle:

            # Results header
            header = "{0},{1},{2},{3},{4}".format("Iteration", "Pitch", "Plunge", "CL", "CDi")
            for y in self.wing.y_cp:
                header += ",{0}".format(y)
            print(header, file=results_file_handle)

            # Loop through iterations
            k = 0
            x_curr = 0.0
            t = 0.0
            while k<self._max_iter:

                # Increment
                k += 1

                # Determine pitch angle
                a = self._a_mean+self._a_A*np.sin(2.0*np.pi*self._w*t+self._phi_a)
                u_n = np.zeros(3)
                u_n[2] = -np.cos(a)
                u_n[0] = -np.sin(a)

                # Determine control and node points
                cp = np.zeros((self.wing.N, 3))
                cp[:,0] = x_curr
                cp[:,1] = self.wing.y_cp
                nodes = np.zeros((self.wing.N+1, 3))
                nodes[:,0] = x_curr
                nodes[:,1] = self.wing.y_node
                if not self.wing.planar:
                    raise IOError

                # Update wake points
                self.wing.wake_points[k,:,:] = np.copy(nodes)

                # Determine plunge velocity
                V_P = 2.0*np.pi*self._z_A*self._w*np.cos(2.0*np.pi*self._w*t+self._phi_z)

                # Determine V_sig due to freestream and plunging
                V_sig[:,0] = -self._V_inf
                V_sig[:,2] = V_P
                V_inf = np.copy(V_sig)
                V_inf_mag = vec_norm(V_inf)

                # Add effect of induced velocity from previous iterations
                v_mji = self.wing.get_influences(cp, k)
                V_sig += np.einsum('ijkl,jk->il', v_mji[:,:k-1,:,:], self.gamma[:k-1,:])

                # Calculate V_sig magnitudes
                V_sig_mag = vec_norm(V_sig)
                V_sig_mag_2 = V_sig_mag**2

                # Calculate b vector
                b = V_inf_mag**2*self.wing.CLa*(vec_inner(V_sig, u_n)/V_inf_mag-self.wing.aL0)*self.wing.dS

                # Assemble A matrix
                A = np.zeros((self.wing.N, self.wing.N))
                A[np.diag_indices(self.wing.N)] = 2.0*vec_norm(vec_cross(V_inf, self.wing.dl))
                A -= V_inf_mag[:,np.newaxis]*self.wing.CLa*self.wing.dS*np.einsum('ikl,l->ik', v_mji[:,k-1,:,:], u_n)

                # Solve for gamma
                self.gamma[k-1,:] = np.linalg.solve(A, b)

                # Plot distribution
                if False:
                    plt.figure()
                    plt.plot(self.wing.y_cp, self.gamma[k-1,:])
                    plt.plot(self.wing.y_cp[::-1], self.gamma[k-1,:])
                    plt.show()

                # Determine force coefficients
                V_sig += np.einsum('ikl,k->il', v_mji[:,k-1,:,:], self.gamma[k-1,:])
                dF = vec_cross(V_sig, self.wing.dl)*self.gamma[k-1,:][:,np.newaxis]
                F = np.sum(dF, axis=0)
                CL = -F[2]/(0.5*self._V_inf**2*self.wing.S)
                # Save CL to array for integration?
                CDi = -F[0]/(0.5*self._V_inf**2*self.wing.S)
                # Save CD to array for intigration?

                # Write results to file
                line = "{0},{1},{2},{3},{4}".format(k, np.degrees(a), V_P, CL, CDi)
                for gamma in self.gamma[k-1]:
                    line += ",{0}".format(gamma)
                print(line, file=results_file_handle)

                # Output
                if verbose: print("{0:<20}{1:<20}{2:<20}{3:<20}{4:<20}".format(k, np.degrees(a), V_P, CL, CDi))

                # Update time and x position
                t += self._dt
                x_curr += self._dt*self._V_inf
            
            # Integrate CL, CD over one period (using file or array stored values)
            '''
            T = 1/self._w
            CL_integ, CD_integ = integrate CL and CD from [self._V_inf*(t-T), self._V_inf*t]
            # period-averaged CL, CD
            CL_T_avg = CL_integ/(self._V_inf*T)
            CL_D_avg = CD_integ/(self._V_inf*T)

            '''




    def export_vtk(self, filename):
        """Writes the case results to a .vtk file.

        Parameters
        ----------
        filename : str
            Name of file to write results to. Must be '.vtk'.
        """

        pass