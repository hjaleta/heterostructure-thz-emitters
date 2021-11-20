import numpy as np
from time import perf_counter
from math import ceil
from EfieldSimulation.help_functions import (get_J, get_J_t, get_r_tuple, 
            back_diff1, back_diff2, interpolation, calc_delta_t)

class Simulation:

    """
        This class contains the entire system. It stores the Mediums, Vacuum, and Efield. It utilizes Jefimenko's equations to model
        the heterostructure and calculate the E-field.
        ...
        ----------
        Attributes
        ----------
            name: str
                name of the simulation. It's convenient to use something that distinguishes
                it from other simulations, like the unique parameters it uses.

            spin_flux: 2d numpy array
                this is the main input to the simulation. It contains the data for the 
                spin current in the mediums. Dimension (time, z)

            sim_p: dict
                Parameters for the simulation. Contains time and spatial discretization, dimensions, 
                specification of which Jefimenko's equations terms that should be used, specifications of 
                region excited by the laser pulse. 
            
            med_p: dict
                Parameters passed to the Medium objecs
            
            vac_p: dict
                Parameters passed to the Vacuum object
            
            time: 1d numpy array
                A vector with the time points
            
            interfaces: list
                z-coordinates of the interfaces between Mediums and Vacuum
            
            mediums: List[Medium]
                The Medium objects of the simulation
            
            vacuum: Vacuum
                The Vacuum object of the simulation

        -------
        Methods
        -------

            The unmentioned methods are used to initialize the class.

            run(print_percent)
                Method that run the entire simulation. If print_percent == True the method
                also prints the progress
            
            step()
                Calculate one timestep of the simulation
            
            save_Ex(path):
                saves the E-field(s) of the simulation as a ndarray txt-file to the specified path
    """

    def __init__(self, spin_flux, sim_params, medium_params, vacuum_params, name = ""):
        self.name = name
        self.spin_flux = spin_flux
        self.sim_p = sim_params
        self.jef_terms = self.sim_p["jef_terms"]
        self.dt = self.sim_p["dt"]
        self.med_p = medium_params
        self.vac_p = vacuum_params
        self.time = self.format_time()
        self.interfaces = self.get_interfaces()
        self.mediums = self.get_mediums()
        self.vacuum = self.get_vacuum()

    def format_time(self):
        time = np.arange(0, self.spin_flux.shape[0], self.dt)
        return time
    
    def get_interfaces(self):
        interfaces = [p["z"][1] for p in self.med_p]
        return interfaces

    def get_mediums(self):
        mediums = []
        
        # Loop over mediums
        for medium in self.med_p: 
            z_0, z_end = medium["z"] # Get medium z-coordinates
            medium_spin_flux = self.spin_flux[:,z_0:z_end] # Extract the spin flux using the coordinates
            mediums.append(Medium(medium, self.sim_p, medium_spin_flux)) # Create Medium
        return mediums

    def get_vacuum(self):
        N_time = len(self.time)
        vac = Vacuum(self.vac_p, self.sim_p, N_time, self.mediums, self.interfaces) # Create Vacuum
        return vac
                            
    def run(self, print_percent = False):
        first_time = perf_counter() #Time the simulation
        
        if print_percent: # If we want progress to be printed
            percent_step = ceil((len(self.time)*print_percent)/100)
            percent = -print_percent
        else:
            percent_step = 0
        max_shift = ceil(self.vacuum.max_delta_t) # The maximum time delay in the simulation
        for t_i, t in enumerate(self.time): # Iterate over each timepoint
            
            self.step(t_i, max_shift) # One timestep
            if print_percent and t_i % percent_step == 0:
                percent += print_percent
                print(f"Simulation {self.name}\n{percent} % completed")
                current_time = perf_counter()
                print(f"Total simulation time: {current_time - first_time:.0f} s\n----------")
        current_time = perf_counter()
        print(f"Simulation {self.name} completed in {current_time-first_time:.0f} s\n----------")

    def step(self, t_i, max_shift):
        dt = self.sim_p["dt"]
        # t = t_i * dt
        dV = self.sim_p["dx"] * self.sim_p["dy"] * self.sim_p["dz"] # The volume elemnt for each discrete cell
        t_0 = t_i - max_shift - 2   # The earliest time we take into consideration
        t_0 = 0 if t_0 < 0 else t_0 # Shift earliest time to 0 if its less than 0
        t_len = t_i - t_0 + 1 # This is how many time points we take inte consideration for one step
        
        J_t_factor = 16022   # prefactor of Jefimenkos equations, converting from particle, nm, fs
                                    # to Coulomb, meter, second

        for m_i, medium in enumerate(self.mediums): # Iterate over Mediums
            Nx, Ny, Nz = medium.N_points


            if medium.params["use_rho"]:
                pass

            if medium.params["use_Jx"]: # This code block handles the ISHE - induced charge current

                # Get the excitation region array for each time point 
                exc_reg = np.tile(medium.exc_region, (t_len, 1, 1, 1)) 

                # Extract the spin current from the medium for the relevant time points
                spin_flux = medium.spin_flux[t_0:t_i+1,:].reshape((t_len, 1, 1, Nz))

                # Calculate the transverse charge flux in each excited region,
                #  with the material's spin hall angle and the spin flux
                charge_flux = exc_reg * spin_flux * medium.params["theta"]

                for E in self.vacuum.E_fields: # Iterate over all E-fields
                    if E.delta_r[m_i] is None: # If Medium does not contribute to E-field, skip it 
                        continue
                    Ex = 0 # Initiate the x-component of the E-field
                    delta_t, delta_r_tot = E.delta_t[m_i], E.delta_r_tot[m_i] # Fetch time delays and source vectors
                    
                    for index, t_shift in np.ndenumerate(delta_t): # Iterate over all discrete cells in the medium
                        x_i, y_i, z_i = index # Get x,y,z indices
                        r = delta_r_tot[index] # Get distance to cell
                        t_i_ret = (t_len - 1 - t_shift)/dt # Calculate the retarded time index for the cell

                        # Get the charge flux for a particular cell, for the considered time points
                        charge_flux_xyz = charge_flux[:,x_i, y_i, z_i].copy().flatten() 

                        # Calculate the derivative of the charge flux
                        J_t = get_J_t(t_i_ret, dt, charge_flux_xyz)
                        
                        # Add the contribution from the cell tot the total E-field
                        Ex += -J_t*J_t_factor*dV/r
                    
                    E.Ex[t_i] = Ex # Add the calculated E-field to the result array
                    

    def save_Ex(self, path): 
        all_arrs = []
        for E in self.vacuum.E_fields:
            E_data = E.Ex.copy().reshape((1,-1))
            all_arrs.append(E_data)
        arr = np.array(all_arrs)
        np.savetxt(path, arr)

class Medium:

    """
        This class represents one component of the heterostructure. 
        ...

        ----------
        Attributes
        ----------
            params: dict:
                Contains the parameters of the medium: material, spatial dimensions, which quantities
                to be used in the simulation, the spin hall angle theta and refractive index n.
            
            spin_flux: 2d numpy array
                the spin flux in this particular medium. Dimensions (time, z)
            
            N_points: tuple: (Nx, Ny, Nz)
                The number of discrete units of the material per respective spatial dimension
            
            exc_region: 2d numpy array
                The region of the Medium that is excited by the laser, and therefore also carries spin_current. 
                The array consists of ones in the excited region, and zeros in the other. Dimension (x, y)

            use: bool
                bool that states whether the medium is used or not in the simulation, i.e. 
                whether it contributes to the Efield

        ------- 
        Methods
        -------
            All methods are used to initialize the class
    """

    def __init__(self, params, sim_params, spin_flux):
        self.params = params
        self.spin_flux = spin_flux
        self.format_domain(sim_params)
        self.N_points = self.get_N_points(sim_params)
        self.exc_region = self.get_excitation_region(sim_params)
        self.use = self.get_usage()
    
    def format_domain(self, sim_params):
        # This function shifts the coordinates of the boundaries of the medium
        # The new coordinates represent the center of each cell, instead of the boundaries of them
        for coord in ["x", "y", "z"]:
            old_min, old_max = self.params[coord]
            shift = sim_params[f"d{coord}"]/2
            self.params[coord] = (old_min + shift, old_max - shift)

    def get_N_points(self, sim_p):
        # This function returns the dimensions of the discretized medium
        N_points = lambda coord: int((self.params[coord][1] - self.params[coord][0])/sim_p[f"d{coord}"] + 1 )
        Nx = N_points("x")
        Ny = N_points("y")
        Nz = N_points("z")
        return (Nx, Ny, Nz)

    def get_excitation_region(self, sim_p):
        # This function returns an array with zeros and ones. The ones represent the exciation region
        # where the laser hits the metal => where we have excied electrons => where we have spin flux
        Nx, Ny, Nz = self.N_points
        exc_region = sim_p["exc_region"]
        x_min, x_max = self.params["x"]
        y_min, y_max = self.params["y"]
        dx = sim_p["dx"]
        dy = sim_p["dy"]
        X = np.arange(x_min, x_max + dx, dx)
        Y = np.arange(y_min, y_max + dy, dy)
        if (len(X) != Nx) or (len(Y) != Ny): 
            raise ValueError("Dimensions of excitation region is inconsistent")

        if exc_region["shape"] == "rectangle": # This code block calculates a rectangular excitation region
            x_side, y_side = exc_region["x_side"]/2, exc_region["y_side"]/2 # x-y coordinates
            l1 = []
            for x_i in X:
                l2 = []
                for y_i in Y:
                    if (abs(x_i) <= x_side) and (abs(y_i) <= y_side): # Check if point is within the rectangle
                        l2.append(1)
                    else:
                        l2.append(0)
                l1.append(l2)

        elif exc_region["shape"] == "circle": # This code block calculates a circular excitation region
            r_squared = exc_region["radius"]**2  # Radius squared
            l1 = []
            for x_i in X:
                l2 = []
                for y_i in Y:
                    if x_i**2 + y_i**2 < r_squared: #Check if point is within circle

                        l2.append(1)
                    else:
                        l2.append(0)
                l1.append(l2)
        else:
            raise ValueError("Invalid shape of excitation region, use 'circle' or 'square'")
        region = np.array([l1]*Nz) # Add z-dimension and create array
        region = np.moveaxis(region, 0, 2) # Fix compatibility of dimension 
        # print(region.shape)
        return region

    def get_usage(self):
        # This function checks if the Medium contributes to the E-field
        for coord in ["x", "y", "z"]:
            if self.params[f"use_J{coord}"]:
                return True
        if self.params["use_rho"] or self.params["use_field"]:
            return True
        else:
            return False

class Vacuum:
    """
        This class represents the vacuum next to the heterostructure. 

        -----------
        Attributes
        -----------
            params: dict:
                Contains the parameters of the vacuum: different z values where the E-field should be calculated,
                and which components of the field that should be calculated
            
            z: list
                list of the 
            
            N_points: tuple: (Nx, Ny, Nz)
                The number of discrete units of the material per respective spatial dimension
            
            exc_region: 2d numpy array
                The region of the Medium that is excited by the laser, and therefore also carries spin_current. 
                The array consists of ones in the excited region, and zeros in the other. Dimension (x, y)

            use: bool
                bool that states whether the medium is used or not in the simulation, i.e. 
                whether it contributes to the Efield
        
        -------
        Methods
        -------
            All methods are used to initialize the class
    """
    def __init__(self, params, sim_params, N_time, mediums, interfaces):
        self.params = params
        self.z = self.get_z_points()
        self.E_fields = self.get_E_fields(sim_params, mediums, N_time, interfaces)
        self.max_delta_t = self.get_max_delta_t()
        
    def get_z_points(self):
        return self.params["z"]

    def get_E_fields(self, sim_params, mediums, N_time, interfaces):
        E_fields = []
        
        for z in self.z: # Iterate over the distances we want to calculate the E-field in
            E = E_field(z, sim_params, mediums, N_time, interfaces) # Create E_field
            E_fields.append(E)
        return E_fields

    def get_max_delta_t(self):
        # This function finds the biggest time retardation for the system
        maxes =[]
        for E in self.E_fields:
            maxes.append(E.delta_t_max)
        return max(maxes)

class E_field:
    """
        This class represents the E-field in one point.

        ----------
        Attributes
        ----------

            z: int or float
                The z-coordinate of the E-field
            
            Ex: 1d numpy array:
                Values of the electric field. Initialized as zeros before the simulation runs.
            
            delta_r: list[4d numpy array or None]:
                This list holds either None objects or 4d arrays. Every element in the list corresponds to a Medium in the
                simulation. If the Medium doesn't have an explicit contribution to the E-field (i.e. contributes to some
                term in the Jefimenko equations), the element is None.
                
                If the Medium on the other hand contributes, the element is a 4d numpy array. These arrays hold the norm of
                every source vector corresponding to the z-coordinate of the E-field. The norm is split into several parts,
                one for each media the source vector passes through. The dimension of the array is: 
                (Nx, Ny, Nz, number of mediums + 1)

                Example: The source vector origins from a point in the 2nd of two media. It does not pass through Medium 1,
                the part of it inside Medium 2 is 2.5 nm long and the part in the Vacuum is 10 000 nm long. The corresponding
                source vector would then hold the values [0, 2.5, 10 000]
            
            delta_r_tot: list[3d numpy array or None]:
                This is the same as the delta_r, but the norm is not split into several parts, but only the total one.

                Dimension of the array: (Nx, Ny, Nz)

                Example: The source vector origins from a point in the 2nd of two media. It does not pass through Medium 1,
                the part of it inside Medium 2 is 2.5 nm long and the part in the Vacuum is 10 000 nm long. The corresponding
                source vector would then be 10 002.5

            delta_t: list[3d numpy array or None:]
                This list holds either None objects or 3d arrays. Every element in the list corresponds to a Medium in the
                simulation. If the Medium doesn't have an explicit contribution to the E-field (i.e. contributes to some
                term in the Jefimenko equations), the element is None.
                
                If the Medium on the other hand contributes, the element is a 3d numpy array. These arrays hold the time
                that it would take for the E-field to traverse from the source point to the place where the E-field 
                is to be calculated. The time corresponds to the source vectors in self.delta_r. The time values are
                calculated with respect to the refractive indices of the material. 
                Dimension of the array is: (Nx, Ny, Nz)

                Example: The source vector origins from a point in the 2nd of two media. It does not pass through Medium 1,
                the part of it inside Medium 2 is 2.5 nm long and the part in the Vacuum is 10 000 nm long. The refractive idex 
                of Medium 2 is 100. Then the time value would be (2.5 * 100 + 10 000)*1/c

        -------   
        Methods
        -------

            All methods are used to initialize the class
    """
    def __init__(self, z, sim_params, mediums, N_time, interfaces):
        self.z = z
        self.Ex = np.zeros(N_time)
        self.delta_r = self.get_delta_r(sim_params, mediums, interfaces)
        self.delta_r_tot = self.get_delta_r_tot()
        self.delta_t, self.delta_t_max = self.get_delta_t(mediums)

    def get_delta_r(self, sim_params, mediums, interfaces):
        # List of array or NoneType correspondning to each medium. If array, it contains source vectors 
        # for the E_field. The source vectors consist of multiple numbers, each corresponding to 
        # the distance in respective Medium / Vacuum
        all_delta_r = []
        for m_i, medium in enumerate(mediums): # Iterate over Mediums
            if medium.use: # If the Medium contributes to the E_field

                # Get spatial properties of Medium
                x_min, x_max = medium.params["x"]
                y_min, y_max = medium.params["y"]
                z_min, z_max = medium.params["z"]
                Nx, Ny, Nz = medium.N_points
                dx, dy, dz =  sim_params["dx"], sim_params["dy"], sim_params["dz"]

                # Allocate space for source vectors
                delta_r = np.zeros((Nx,Ny,Nz, len(mediums)+1))

                for (x_i,y_i,z_i), val in np.ndenumerate(delta_r[:,:,:,0]): # Iterate over all cells in Medium

                    # Calculate source vector coordinates
                    x_s = x_min + x_i * dx
                    y_s = y_min + y_i * dy
                    z_s = z_min + z_i * dz

                    # Calculate distances in source vector. 
                    # Every distance corresponds to one Medium or Vacuum
                    r = get_r_tuple(x_s, y_s, z_s, self.z, interfaces)

                    # Put source vector into array
                    delta_r[x_i,y_i,z_i,:] = r
                all_delta_r.append(delta_r)
            else: # If no contribution from Medium, just append None
                all_delta_r.append(None)
        return all_delta_r

    def get_delta_r_tot(self):
        # Calculate the total distance of the source vectors
        all_delta_r_tot = []
        for delta_r in self.delta_r: # Iterate over Mediums
            if delta_r is None:
                all_delta_r_tot.append(None)
            else:
                # Allocate space for array
                r_array = np.zeros(delta_r[:,:,:,0].shape)
                # Iterae over source vectors
                for (xi,yi,zi), _ in np.ndenumerate(delta_r[:,:,:,0]):
                    r_tot = np.sum(delta_r[xi,yi,zi,:]) # Sum contributions to total source vector
                    r_array[xi,yi,zi] = r_tot
                all_delta_r_tot.append(r_array)
        return all_delta_r_tot

    def get_delta_t(self, mediums):
        n = []
        for medium in mediums: # Fetch refravctive indices
            n.append(medium.params["n"])
        n.append(1) # Append the Vacuum refractive index
        all_delta_t = []
        delta_t_maxes = []
        for arr in self.delta_r: # Iterate over source vectors
            if arr is None: # If no contribution
                all_delta_t.append(None)
            else:
                Nx,Ny,Nz,Nm = arr.shape 
                delta_t = np.zeros((Nx,Ny,Nz)) # Allocate space for time delays
                for (x_i, y_i, z_i), _ in np.ndenumerate(delta_t): # Iterate over source vectors
                    tup = arr[x_i,y_i,z_i,:].copy() # Fetch source vector
                    t = calc_delta_t(tup, n) # Calculate time retardation
                    delta_t[x_i,y_i,z_i] = t # Put time retardation into array
                t_min = np.amin(delta_t) # Calculate the smallest time delay
                delta_t -= np.ones(delta_t.shape)*t_min # Shift time delays so that the smallest time delay is zero
                all_delta_t.append(delta_t)
                delta_t_maxes.append(np.amax(delta_t))
        delta_t_max = max(delta_t_maxes)

        return all_delta_t, delta_t_max

if __name__ == "__main__":
    pass
    