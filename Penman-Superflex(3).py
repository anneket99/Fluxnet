#!/usr/bin/env python
# coding: utf-8

# In[76]:


from superflexpy.implementation.elements.hbv import PowerReservoir
from superflexpy.implementation.elements.gr4j import UnitHydrograph1
from superflexpy.implementation.root_finders.pegasus import PegasusPython
from superflexpy.framework.element import ODEsElement
from superflexpy.implementation.numerical_approximators.implicit_euler import ImplicitEulerPython
from superflexpy.framework.unit import Unit
from superflexpy.framework.node import Node
from superflexpy.implementation.elements.structure_elements import Junction, Splitter, Transparent
import numba as nb


# ![grafik.png](attachment:grafik.png)

# In[77]:


solver_python = PegasusPython()
approximator = ImplicitEulerPython(root_finder=solver_python)


# In[78]:


class RootZoneStorage(ODEsElement):
    
    #Input: Precipitation, Evapotranspiration
    #Output: q_ex
    def __init__(self, parameters, states, approximation, id):

        ODEsElement.__init__(self,
                             parameters=parameters,
                             states=states,
                             approximation=approximation,
                             id=id)

        self._fluxes_python = [self._fluxes_function_python]

        if approximation.architecture == 'python':
            self._fluxes = [self._fluxes_function_python]

    def set_input(self, input):
        self.input = {'P': input[0],
                     'PET': input[1]}   #Precipitation
                     
    def get_output(self, solve=True):

        if solve:
            self._solver_states = [self._states[self._prefix_states + 'S0']]

            self._solve_differential_equation()
            # Update the state
            self.set_states({self._prefix_states + 'S0': self.state_array[-1, 0]})

        fluxes = self._num_app.get_fluxes(fluxes=self._fluxes_python,
                                          S=self.state_array,
                                          S0=self._solver_states,
                                          dt=self._dt,
                                          **self.input,
                                          **{k[len(self._prefix_parameters):]: self._parameters[k] for k in self._parameters},
                                          )
        return [-fluxes[0][2]]

    def get_AET(self):
        ''' Calculates Actual Evapotranspiration'''

        try:
            S = self.state_array
        except AttributeError:
            message = '{}get_aet method has to be run after running '.format(self._error_message)
            message += 'the model using the method get_output'
            raise AttributeError(message)

        fluxes = self._num_app.get_fluxes(fluxes=self._fluxes_python,
                                          S=S,
                                          S0=self._solver_states,
                                          dt=self._dt,
                                          **self.input,
                                          **{k[len(self._prefix_parameters):]: self._parameters[k] for k in self._parameters},
                                          )
        return [- fluxes[0][1]]


    @staticmethod
    def _fluxes_function_python(S, S0, ind, P, PET,S_max,dt):
        ''' Differential Equation to calculate S(Root Zone Storage) from P and PET'''
        if ind is None:
            
            # Calculate Potential Evapotranpiration
            if(PET[ind] >= 0):
            
                E_a = PET[ind]
            else:
                E_a = 0
            # Calculate excess
            
            if (S + P >= S_max):
        
                q_ex = S + P -S_max
            else:
                q_ex = 0

        
            return (P-E_a-q_ex)*dt+S_0
           
        else:
            # Calculate Potential Evapotranpiration
            if(PET[ind] >= 0):
                E_a = PET[ind]
            else:
                E_a = 0
                
             # Calculate excess
            if (S[ind] + P[ind] >= S_max[ind]):
                q_ex = S[ind] + P[ind]-S_max[ind]
            else:
                q_ex = 0

    return (P[ind]-E_a-q_ex)*dt+S_0


# In[80]:


class MoistureDeficit(ODEsElement):
    'Calculates Moisture Deficit from Excess'

    def __init__(self, parameters, states, approximation, id):

        ODEsElement.__init__(self,
                             parameters=parameters,
                             states=states,
                             approximation=approximation,
                             id=id)

        self._fluxes_python = [self._fluxes_function_python]

        if approximation.architecture == 'python':
            self._fluxes = [self._fluxes_function_python]

    def set_input(self, input):

        self.input = {'P': input[0]}
                     

    def get_output(self, solve=True):

        if solve:
            self._solver_states = [self._states[self._prefix_states + 'S_def']]

            self._solve_differential_equation()

            # Update the state
            self.set_states({self._prefix_states + 'S_def_0': self.state_array[-1, 0]})

        fluxes = self._num_app.get_fluxes(fluxes=self._fluxes_python,
                                          S_def=self.state_array,
                                          S_def_0=self._solver_states,
                                          dt=self._dt,
                                          **self.input,
                                          **{k[len(self._prefix_parameters):]: self._parameters[k] for k in self._parameters},
                                          )
        return [-fluxes[0][1]]

   
    @staticmethod
    def _fluxes_function_python(S_def, S_def_0,ind,gamma,PET,dt):
        if ind is None:
            
            # Calcualte soil moisture excess (u2)
            if(S_def <= 0):
                u2 = q12

            if(S_def - q12 <= 0):
                u2 = q12-S_def

            else:
                u2 = 0

            S_def = (E_t + u2 - q12)*dt + S_def_0
            # Calculate Evapotranspiration
            if(S_rz <= 0):
                E_t = gamma*PET
            else: 
                E_t = 0
                return E_t

            S_def = (E_t + u2 - q12)*dt + S_def_0

            return S_def

        if ind:
            if(S_def <= 0):
                u2 = q12[ind]

            if(S_def[ind] - q12 <= 0):
                    u2 = q12[ind]-S_def[ind]

            else:
                u2 = 0

            S_def[ind] = (E_t[ind] + u2 - q12[ind])*dt + S_def_0

            if(S_rz[ind] <= 0):
                E_t = gamma*PET[ind]
            else: 
                E_t = 0
                
                return E_t

            S_def = (E_t + u2 - q12)*dt + S_def_0
            return S_def
            
            


# In[81]:


##Current Storage in the Root Zone
class RootingReservoir(ODEsElement):
    ''' Calculate Rooting Reservoir from u1 and u2'''

    def __init__(self, parameters, states, approximation, id):

        ODEsElement.__init__(self,
                             parameters=parameters,
                             states=states,
                             approximation=approximation,
                             id=id)

        self._fluxes_python = [self._fluxes_function_python]

        if approximation.architecture == 'numba':
            self._fluxes = [self._fluxes_function_numba]
        elif approximation.architecture == 'python':
            self._fluxes = [self._fluxes_function_python]

    def set_input(self, input):

        self.input = {'P': input[0],
                      'PET': input[1]}

    def get_output(self, solve=True):

        if solve:
            self._solver_states = [self._states[self._prefix_states + 'C_res_0']]

            self._solve_differential_equation()

            # Update the state
            self.set_states({self._prefix_states + 'C_res_0': self.state_array[-1, 0]})

        fluxes = self._num_app.get_fluxes(fluxes=self._fluxes_python,
                                          C_res=self.state_array,
                                          C_res_0=self._solver_states,
                                          dt=self._dt,
                                          **self.input,
                                          **{k[len(self._prefix_parameters):]: self._parameters[k] for k in self._parameters},
                                         )
        return [-fluxes[0][2]]


    # PROTECTED METHODS

    @staticmethod
    def _fluxes_function_python(C_res, C_res_0, ind,S,P,S_max,dt):
        if ind is None: 
            # Check if Rooting Reservoir is at maximum
            if (S + P >= S_max):
        
                q_ex = S + P -S_max
            else:
                q_ex = 0
                
            # Calculate Outflow
            Q = k1*C_res
            C_res = C_res_0 + (q_ex -Q)*dt
    
            return C_res
        else:
            # Check if Rooting Reservoir is at maximum
            if (S[ind] + P[ind]>= S_max):
        
                q_ex = S + P -S_max
            else:
                q_ex = 0
                
            # Calculate Outflow 
                Q = k1*C_res[ind]
                C_res = C_res_0 + (q_ex-Q)*dt
            
            return C_res


# In[ ]:





# In[82]:


root_finder = PegasusPython()  # Use the default parameters
numerical_approximation = ImplicitEulerPython(root_finder)

# Set up parameters
rootzone_storage = RootZoneStorage(parameters={'Smax': 200.0,},
                       states={'S0': 50.0},
                       approximation=numerical_approximation,
                       id='rs')

moisture_deficit = MoistureDeficit(parameters={ 'gamma' : 0.08},
                                    states={'S_def_0': 80.0},
                                    approximation=numerical_approximation,
                                    id='md')
rooting_reservoir = RootZoneStorage(parameters={'Smax': 200.0},
                                   states = {'C_res_0': 100.0},
                                   approximation=numerical_approximation,
                                    id='rr')
# Splitter to split qex into u1 and q12
s = Splitter(
    weight=[[1-0.15], [0.15]],
    direction=[[0], [0]],
    id='S'
)
# Junction combines u1 and u2
j = Junction(
    direction=[[0,0],
               [0,0]],
    id='J'
)

t = Transparent(id= 'T')

# Set up model (combine layers)
model = Unit(layers=[
    [rootzone_storage],
    [s],
    [moisture_deficit, t],
    [j],
    [rooting_reservoir]
],
             id = 'model')


# In[83]:


import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})

# create Precipitation Data

# Fix the seed
SEED = 2
rng = np.random.RandomState(seed=SEED)
#Generate the input
P = np.zeros(100)
P[:10] = rng.randint(10, size=10)
P[25:30] = rng.randint(20, size=5)
P[40:60] = rng.randint(5, size=20)
P[80:83] = rng.randint(30, 50, size=3)
E = np.ones_like(P) * 2.0 # We use a constant PET
# Assign the input
model.set_input([P, E])
# Set the timestep
model.set_timestep(1.0)
# Run the model
model.reset_states()
output = model.get_output()

# Inspect internals
rs_out = model.call_internal(id='rs', method='get_output', solve=False)[0]
md_out = model.call_internal(id='md', method='get_output', solve=False)[0]
rr_out = model.call_internal(id='rr', method='get_output', solve=False)[0]
rs_e = model.call_internal(id='rs', method='get_AET')[0]
rs_s = model.get_internal(id='rs', attribute='state_array')[:, 0]
md_s = model.get_internal(id='md', attribute='state_array')[:, 0]
rr_s = model.get_internal(id='rr', attribute='state_array')[:, 0]

# Plot
fig, ax = plt.subplots(3, 1, figsize=(20, 12), sharex=True)
ax[0].bar(x=np.arange(len(P)), height=P, color='royalblue', label='P')
ax[0].plot(np.arange(len(P)), E, lw=2, color='gold', label='PET')
ax[0].legend()
ax[0].set_ylabel('Inputs [mm/day]')
ax[0].grid(True)
ax[1].plot(np.arange(len(P)), output[0], lw=2, label='Total outflow')
ax[1].plot(np.arange(len(P)), rs_e, lw=2, label='AET')
ax[1].plot(np.arange(len(P)), rs_out, lw=2, label='Outflow Root Zone Storage')
ax[1].plot(np.arange(len(P)), md_out, lw=2, label='Outflow Moisture Deficit')
ax[1].plot(np.arange(len(P)), rr_out, lw=2, label='Outflow ')
ax[1].set_xlabel('Time [days]')
ax[1].set_ylabel('Flows [mm/day]')
ax[1].legend()
ax[1].grid(True)
ax[2].plot(np.arange(len(P)), rs_s, lw=2, label='State Root Zone Storage')
ax[2].plot(np.arange(len(P)), md_s, lw=2, label='State Moisture Deficit')
ax[2].plot(np.arange(len(P)), rr_s, lw=2, label='State Rooting Reservoir')
ax[2].set_xlabel('Time [days]')
ax[2].set_ylabel('Storages [mm]')
ax[2].legend()
ax[2].grid(True)
pass


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




