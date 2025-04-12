'''
Advances.py :

Under Developing







'''

#" IN GOD WE TRUST, ALL OTHERS MUST BRING DATA"
#                                               -W. Edwards Deming
#------------------------------------------------------------------------------
# Copyright 2023 The Gamlab Authors. All Rights Reserved.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#------------------------------------------------------------------------------
''' 
The Scientific experimental simulation library 
-------------------------------------------------------------------------------
Graphen & Advanced Material Laboratory 

it aimes to provide new scientist to use data,simlation, prepared data 
and Artificial intelligence models.

See http://gamlab.aut.ac.ir for complete documentation.
'''
__doc__='''

@author: Ali Pilehvar Meibody (Alipilehvar1999@gmail.com)

                                         888                    888
 .d8888b    .d88b.     88888b.d88b.      888         .d88b.     888
d88P"      d88""88b    888 "888 "88b     888        d88""88b    88888PP
888  8888  888  888    888  888  888     888        888  888    888  888
Y88b.  88  Y88..88PP.  888  888  888     888......  Y88..88PP.  888  888
 "Y8888P8   "Y88P8888  888  888  888     888888888   "Y88P8888  88888888  


@Director of Gamlab: Professor M. Naderi (Mnaderi@aut.ac.ir)    

@Graphene Advanced Material Laboratory: https://www.GamLab.Aut.ac.ir


@Co-authors: 
'''




'''
#import-----------------------------------------
import math
import statistics
import cmath
import random
import numpy as np
import pandas as pd
import matplotlib as plt





from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

class NavierStokesSolver:
    def __init__(self, domain_size, viscosity, density):
        self.domain_size = domain_size
        self.viscosity = viscosity
        self.density = density

        # Define grid
        self.grid_size = (domain_size, domain_size)
        self.velocity_field = np.zeros((2, *self.grid_size))
        self.pressure_field = np.zeros(self.grid_size)

    def solve_step(self, dt):
        # Implement a simple explicit finite difference method for illustration
        # Update velocity field
        self.velocity_field += dt * self.calculate_velocity_change()

        # Update pressure field
        self.pressure_field = self.calculate_pressure()

    def calculate_velocity_change(self):
        # Placeholder function for calculating velocity change
        # Implement based on simplified Navier-Stokes equations
        # This is a basic example, and you may need more advanced methods for accuracy
        velocity_change = np.zeros_like(self.velocity_field)
        return velocity_change

    def calculate_pressure(self):
        # Placeholder function for calculating pressure field
        # Implement based on simplified Navier-Stokes equations
        # This is a basic example, and you may need more advanced methods for accuracy
        pressure = np.zeros_like(self.pressure_field)
        return pressure

    def plot_results(self):
        x, y = np.meshgrid(range(self.domain_size), range(self.domain_size))

        # Plot velocity field
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.quiver(x, y, self.velocity_field[0], self.velocity_field[1], scale=20)
        plt.title('Velocity Field')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')

        # Plot pressure field
        plt.subplot(1, 2, 2)
        plt.imshow(self.pressure_field, cmap='viridis')
        plt.colorbar()
        plt.title('Pressure Field')
        plt.show()
        
    def animate(self, num_frames, dt):
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        def update(frame):
            self.solve_step(dt)
            ax[0].clear()
            ax[0].quiver(x, y, self.velocity_field[0], self.velocity_field[1], scale=20)
            ax[0].set_title('Velocity Field')
            ax[0].set_xlabel('X-axis')
            ax[0].set_ylabel('Y-axis')

            ax[1].clear()
            ax[1].imshow(self.pressure_field, cmap='viridis')
            ax[1].set_title('Pressure Field')

        x, y = np.meshgrid(range(self.domain_size), range(self.domain_size))
        ani = FuncAnimation(fig, update, frames=num_frames, repeat=False)
        plt.show()


'''





#**
'''
# Example Usage:
domain_size = 50
viscosity = 0.01
density = 1.0
dt = 0.01
solver = NavierStokesSolver(domain_size, viscosity, density)

# Time-stepping loop
for _ in range(100):
    solver.solve_step(dt)

# Plot the results
solver.plot_results()

'''




'''
domain_size = 50
viscosity = 0.01
density = 1.0
dt = 0.01
solver = NavierStokesSolver(domain_size, viscosity, density)

# Animate the simulation for 100 frames with a time step of 0.01
solver.animate(num_frames=100, dt=dt)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class NavierStokesSolver:
    def __init__(self, domain_size, viscosity, density):
        self.domain_size = domain_size
        self.viscosity = viscosity
        self.density = density

        # Define grid
        self.grid_size = (domain_size, domain_size)
        self.velocity_field = np.zeros((2, *self.grid_size))
        self.pressure_field = np.zeros(self.grid_size)

        # Initialize with some initial conditions (you may need to modify this based on your problem)
        self.velocity_field[0, :, :] = 1.0  # Initial x-velocity
        self.velocity_field[1, :, :] = 0.0  # Initial y-velocity
    def solve_step(self, dt):
        # Placeholder function for solving the Navier-Stokes equations numerically
        # Implement based on simplified Navier-Stokes equations
        # This is a basic example, and you may need more advanced methods for accuracy

        # Calculate velocity gradients
        du_dx, du_dy = np.gradient(self.velocity_field[0], axis=(1, 2))
        dv_dx, dv_dy = np.gradient(self.velocity_field[1], axis=(1, 2))

        # Calculate Laplacian of velocity field
        d2u_dx2, d2u_dy2 = np.gradient(du_dx, axis=(1, 2)), np.gradient(du_dy, axis=(1, 2))
        d2v_dx2, d2v_dy2 = np.gradient(dv_dx, axis=(1, 2)), np.gradient(dv_dy, axis=(1, 2))

        # Calculate pressure change
        pressure_change = np.zeros_like(self.pressure_field)
        pressure_change[1:-1, 1:-1] = -density * (
            du_dx[1:-1, 1:-1] * (du_dx[1:-1, 1:-1] + dv_dy[1:-1, 1:-1]) +
            du_dy[1:-1, 1:-1] * (dv_dx[1:-1, 1:-1] + dv_dy[1:-1, 1:-1])
        )

        # Update velocity field based on pressure change
        self.velocity_field[0] += dt * (1 / density) * du_dx + dt * viscosity * (d2u_dx2 + d2u_dy2)
        self.velocity_field[1] += dt * (1 / density) * du_dy + dt * viscosity * (d2v_dx2 + d2v_dy2)

        # Update pressure field
        self.pressure_field += dt * pressure_change

        # Enforce boundary conditions (you may need to modify this based on your problem)
        self.velocity_field[:, 0, :] = 0  # Velocity is zero at the bottom boundary
        self.velocity_field[:, -1, :] = 0  # Velocity is zero at the top boundary
        self.velocity_field[:, :, 0] = 0  # Velocity is zero at the left boundary
        self.velocity_field[:, :, -1] = 0  # Velocity is zero at the right boundary

        self.pressure_field[:, 0] = self.pressure_field[:, 1]  # Pressure is constant at the bottom boundary
        self.pressure_field[:, -1] = self.pressure_field[:, -2]  # Pressure is constant at the top boundary
        self.pressure_field[0, :] = self.pressure_field[1, :]  # Pressure is constant at the left boundary
        self.pressure_field[-1, :] = self.pressure_field[-2, :] 
        
    def solve_step2(self, dt):
        # Implement a simple explicit finite difference method for illustration
        # Update velocity field
        self.velocity_field += dt * self.calculate_velocity_change()

        # Update pressure field
        self.pressure_field = self.calculate_pressure()

    def calculate_velocity_change(self):
        # Placeholder function for calculating velocity change
        # Implement based on simplified Navier-Stokes equations
        # This is a basic example, and you may need more advanced methods for accuracy
        velocity_change = np.zeros_like(self.velocity_field)
        return velocity_change

    def calculate_pressure(self):
        # Placeholder function for calculating pressure field
        # Implement based on simplified Navier-Stokes equations
        # This is a basic example, and you may need more advanced methods for accuracy
        pressure = np.zeros_like(self.pressure_field)
        return pressure

    def plot_results(self):
        x, y = np.meshgrid(range(self.domain_size), range(self.domain_size))

        # Plot velocity field
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.quiver(x, y, self.velocity_field[0], self.velocity_field[1], scale=20)
        plt.title('Velocity Field')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')

        # Plot pressure field
        plt.subplot(1, 2, 2)
        plt.imshow(self.pressure_field, cmap='viridis')
        plt.colorbar()
        plt.title('Pressure Field')

    def animate(self, num_frames, dt):
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        def update(frame):
            self.solve_step(dt)
            ax[0].clear()
            ax[0].quiver(x, y, self.velocity_field[0], self.velocity_field[1], scale=20)
            ax[0].set_title('Velocity Field')
            ax[0].set_xlabel('X-axis')
            ax[0].set_ylabel('Y-axis')

            ax[1].clear()
            ax[1].imshow(self.pressure_field, cmap='viridis')
            ax[1].set_title('Pressure Field')

        x, y = np.meshgrid(range(self.domain_size), range(self.domain_size))
        ani = FuncAnimation(fig, update, frames=num_frames, repeat=False)
        plt.show()

# Example Usage:
domain_size = 50
viscosity = 0.01
density = 1.0
dt = 0.01
solver = NavierStokesSolver(domain_size, viscosity, density)
solver.plot_results()
# Animate the simulation for 100 frames with a time step of 0.01
solver.animate(num_frames=100, dt=dt)



import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class NavierStokesSolver:
    def __init__(self, domain_size, viscosity, density):
        self.domain_size = domain_size
        self.viscosity = viscosity
        self.density = density

        # Define grid
        self.grid_size = (domain_size, domain_size)
        self.velocity_field = np.zeros((2, *self.grid_size))
        self.pressure_field = np.zeros(self.grid_size)

        # Initialize with some initial conditions (you may need to modify this based on your problem)
        self.velocity_field[0, :, :] = 1.0  # Initial x-velocity
        self.velocity_field[1, :, :] = 0.0  # Initial y-velocity

    def solve_step(self, dt):
        # Calculate velocity gradients
        du_dx, du_dy = np.gradient(self.velocity_field[0], axis=(1, 2))
        dv_dx, dv_dy = np.gradient(self.velocity_field[1], axis=(1, 2))

        # Calculate Laplacian of velocity field
        d2u_dx2, d2u_dy2 = np.gradient(du_dx, axis=(1, 2)), np.gradient(du_dy, axis=(1, 2))
        d2v_dx2, d2v_dy2 = np.gradient(dv_dx, axis=(1, 2)), np.gradient(dv_dy, axis=(1, 2))

        # Calculate pressure change
        pressure_change = np.zeros_like(self.pressure_field)
        pressure_change[1:-1, 1:-1] = -self.density * (
            du_dx[1:-1, 1:-1] * (du_dx[1:-1, 1:-1] + dv_dy[1:-1, 1:-1]) +
            du_dy[1:-1, 1:-1] * (dv_dx[1:-1, 1:-1] + dv_dy[1:-1, 1:-1])
        )

        # Update velocity field based on pressure change
        self.velocity_field[0] += dt * (1 / self.density) * du_dx + dt * self.viscosity * (d2u_dx2 + d2u_dy2)
        self.velocity_field[1] += dt * (1 / self.density) * du_dy + dt * self.viscosity * (d2v_dx2 + d2v_dy2)

        # Update pressure field
        self.pressure_field += dt * pressure_change

        # Enforce boundary conditions (you may need to modify this based on your problem)
        self.velocity_field[:, 0, :] = 0  # Velocity is zero at the bottom boundary
        self.velocity_field[:, -1, :] = 0  # Velocity is zero at the top boundary
        self.velocity_field[:, :, 0] = 0  # Velocity is zero at the left boundary
        self.velocity_field[:, :, -1] = 0  # Velocity is zero at the right boundary

        self.pressure_field[:, 0] = self.pressure_field[:, 1]  # Pressure is constant at the bottom boundary
        self.pressure_field[:, -1] = self.pressure_field[:, -2]  # Pressure is constant at the top boundary
        self.pressure_field[0, :] = self.pressure_field[1, :]  # Pressure is constant at the left boundary
        self.pressure_field[-1, :] = self.pressure_field[-2, :]  # Pressure is constant at the right boundary

    def initialize_animation(self):
        # Create a figure and axis for the animation
        self.fig, self.ax = plt.subplots(1, 2, figsize=(12, 5))

        # Initialize quiver plot for velocity field
        self.quiver = self.ax[0].quiver([], [], [], [], scale=20)
        self.ax[0].set_title('Velocity Field')
        self.ax[0].set_xlabel('X-axis')
        self.ax[0].set_ylabel('Y-axis')

        # Initialize image plot for pressure field
        self.image = self.ax[1].imshow(self.pressure_field, cmap='viridis')
        self.ax[1].set_title('Pressure Field')

    def update_animation(self, frame):
        # Solve one time step
        self.solve_step(dt)

        # Update quiver plot for velocity field
        self.quiver.set_UVC(self.velocity_field[0], self.velocity_field[1])

        # Update image plot for pressure field
        self.image.set_array(self.pressure_field)

    def animate(self, num_frames, dt):
        # Initialize the animation
        self.initialize_animation()

        # Create the animation
        animation = FuncAnimation(self.fig, self.update_animation, frames=num_frames, interval=50, repeat=False)

        # Show the animation
        plt.show()
        animation.save('navierstoke_test.mp4', writer='ffmpeg', fps=30)

# Example Usage:
domain_size = 50
viscosity = 1
density = 2
dt = 1
solver = NavierStokesSolver(domain_size, viscosity, density)

# Animate the simulation for 100 frames with a time step of 0.01
solver.animate(num_frames=100, dt=dt)


class NavierStokesSolver:
    def __init__(self, domain_size, viscosity, density):
        self.domain_size = domain_size
        self.viscosity = viscosity
        self.density = density

        # Define grid
        self.grid_size = (domain_size, domain_size)
        self.velocity_field = np.zeros((2, *self.grid_size))
        self.pressure_field = np.zeros(self.grid_size)

        # Initialize with some initial conditions (you may need to modify this based on your problem)
        self.velocity_field[0, :, :] = 1.0  # Initial x-velocity
        self.velocity_field[1, :, :] = 0.0  # Initial y-velocity

    def solve_step(self, dt):
        # Calculate velocity gradients
        du_dx, du_dy = np.gradient(self.velocity_field[0], axis=(1, 2))
        dv_dx, dv_dy = np.gradient(self.velocity_field[1], axis=(1, 2))

        # Calculate Laplacian of velocity field
        d2u_dx2, d2u_dy2 = np.gradient(du_dx, axis=(1, 2)), np.gradient(du_dy, axis=(1, 2))
        d2v_dx2, d2v_dy2 = np.gradient(dv_dx, axis=(1, 2)), np.gradient(dv_dy, axis=(1, 2))

        # Calculate pressure change
        pressure_change = np.zeros_like(self.pressure_field)
        pressure_change[1:-1, 1:-1] = -self.density * (
            du_dx[1:-1, 1:-1] * (du_dx[1:-1, 1:-1] + dv_dy[1:-1, 1:-1]) +
            du_dy[1:-1, 1:-1] * (dv_dx[1:-1, 1:-1] + dv_dy[1:-1, 1:-1])
        )

        # Update velocity field based on pressure change
        self.velocity_field[0] += dt * (1 / self.density) * du_dx + dt * self.viscosity * (d2u_dx2 + d2u_dy2)
        self.velocity_field[1] += dt * (1 / self.density) * du_dy + dt * self.viscosity * (d2v_dx2 + d2v_dy2)

        # Update pressure field
        self.pressure_field += dt * pressure_change

        # Enforce boundary conditions (you may need to modify this based on your problem)
        self.velocity_field[:, 0, :] = 0  # Velocity is zero at the bottom boundary
        self.velocity_field[:, -1, :] = 0  # Velocity is zero at the top boundary
        self.velocity_field[:, :, 0] = 0  # Velocity is zero at the left boundary
        self.velocity_field[:, :, -1] = 0  # Velocity is zero at the right boundary

        self.pressure_field[:, 0] = self.pressure_field[:, 1]  # Pressure is constant at the bottom boundary
        self.pressure_field[:, -1] = self.pressure_field[:, -2]  # Pressure is constant at the top boundary
        self.pressure_field[0, :] = self.pressure_field[1, :]  # Pressure is constant at the left boundary
        self.pressure_field[-1, :] = self.pressure_field[-2, :]  # Pressure is constant at the right boundary


    def initialize_animation(self):
        # Create a figure and axis for the animation
        self.fig, self.ax = plt.subplots(1, 2, figsize=(12, 5))

        # Initialize quiver plot for velocity field
        self.quiver = self.ax[0].quiver([], [], [], [], scale=20)
        self.ax[0].set_title('Velocity Field')
        self.ax[0].set_xlabel('X-axis')
        self.ax[0].set_ylabel('Y-axis')

        # Initialize image plot for pressure field
        self.image = self.ax[1].imshow(self.pressure_field, cmap='viridis')
        self.ax[1].set_title('Pressure Field')

    def update_animation(self, frame):
        # Solve one time step
        self.solve_step(dt)

        # Update quiver plot for velocity field
        self.quiver.set_UVC(self.velocity_field[0], self.velocity_field[1])

        # Update image plot for pressure field
        self.image.set_array(self.pressure_field)

    def animate(self, num_frames, dt):
        # Initialize the animation
        self.initialize_animation()

        # Create the animation
        animation = FuncAnimation(self.fig, self.update_animation, frames=num_frames, interval=50, repeat=False)

        # Show the animation

        plt.show()
        animation.save('navierstoke_test.mp4', writer='ffmpeg', fps=30)

# Example Usage:
domain_size = 50
viscosity = 0.01
density = 1.0
dt = 0.01
solver = NavierStokesSolver(domain_size, viscosity, density)

# Animate the simulation for 100 frames with a time step of 0.01
solver.animate(num_frames=100, dt=dt)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Create a simple animation of a sine wave
fig, ax = plt.subplots()
x = np.linspace(0, 2 * np.pi, 100)
line, = ax.plot(x, np.sin(x))

def update(frame):
    line.set_ydata(np.sin(x + frame / 10))  # Update the sine wave
    return line,

ani = FuncAnimation(fig, update, frames=100, interval=50, blit=True)

# Show the animation
plt.show()
ani.save('animation_test.mp4', writer='ffmpeg', fps=30)





import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class NavierStokesSolver:
    def __init__(self, domain_size, viscosity, density):
        self.domain_size = domain_size
        self.viscosity = viscosity
        self.density = density

        # Define grid
        self.grid_size = (domain_size, domain_size)
        self.velocity_field = np.zeros((2, *self.grid_size))
        self.pressure_field = np.zeros(self.grid_size)

        # Initialize with some initial conditions
        self.velocity_field[0, :, :] = 1.0  # Initial x-velocity
        self.velocity_field[1, :, :] = 0.0  # Initial y-velocity

    def solve_step(self, dt):
        # Calculate velocity gradients
        du_dx, du_dy = np.gradient(self.velocity_field[0], axis=(0, 1))
        dv_dx, dv_dy = np.gradient(self.velocity_field[1], axis=(0, 1))

        # Calculate Laplacian of velocity field
        d2u_dx2 = np.gradient(du_dx, axis=(0, 1))
        d2u_dy2 = np.gradient(du_dy, axis=(0, 1))
        d2v_dx2 = np.gradient(dv_dx, axis=(0, 1))
        d2v_dy2 = np.gradient(dv_dy, axis=(0, 1))

        # Calculate pressure change
        pressure_change = np.zeros_like(self.pressure_field)
        pressure_change[1:-1, 1:-1] = -self.density * (
            du_dx[1:-1, 1:-1] * (du_dx[1:-1, 1:-1] + dv_dy[1:-1, 1:-1]) +
            du_dy[1:-1, 1:-1] * (dv_dx[1:-1, 1:-1] + dv_dy[1:-1, 1:-1])
        )

        # Update velocity field based on pressure change
        self.velocity_field[0, 1:-1, 1:-1] += dt * (1 / self.density) * du_dx[1:-1, 1:-1] + dt * self.viscosity * (d2u_dx2[1:-1, 1:-1] + d2u_dy2[1:-1, 1:-1])
        self.velocity_field[1, 1:-1, 1:-1] += dt * (1 / self.density) * du_dy[1:-1, 1:-1] + dt * self.viscosity * (d2v_dx2[1:-1, 1:-1] + d2v_dy2[1:-1, 1:-1])

        # Update pressure field
        self.pressure_field[1:-1, 1:-1] += dt * pressure_change[1:-1, 1:-1]

        # Enforce boundary conditions
        self.velocity_field[:, 0, :] = 0  # Velocity is zero at the bottom boundary
        self.velocity_field[:, -1, :] = 0  # Velocity is zero at the top boundary
        self.velocity_field[:, :, 0] = 0  # Velocity is zero at the left boundary
        self.velocity_field[:, :, -1] = 0  # Velocity is zero at the right boundary

        self.pressure_field[:, 0] = self.pressure_field[:, 1]  # Pressure is constant at the bottom boundary
        self.pressure_field[:, -1] = self.pressure_field[:, -2]  # Pressure is constant at the top boundary
        self.pressure_field[0, :] = self.pressure_field[1, :]  # Pressure is constant at the left boundary
        self.pressure_field[-1, :] = self.pressure_field[-2, :]  # Pressure is constant at the right boundary

    def initialize_animation(self):
        self.fig, self.ax = plt.subplots(1, 2, figsize=(12, 5))

        # Initialize quiver plot for velocity field
        self.quiver = self.ax[0].quiver([], [], [], [], scale=20)
        self.ax[0].set_title('Velocity Field')
        self.ax[0].set_xlabel('X-axis')
        self.ax[0].set_ylabel('Y-axis')

        # Initialize image plot for pressure field
        self.image = self.ax[1].imshow(self.pressure_field, cmap='viridis')
        self.ax[1].set_title('Pressure Field')

    def update_animation(self, frame):
        self.solve_step(dt)

        # Update quiver plot for velocity field
        self.quiver.set_UVC(self.velocity_field[0, 1:-1, 1:-1], self.velocity_field[1, 1:-1, 1:-1])

        # Update image plot for pressure field
        self.image.set_array(self.pressure_field)

    def animate(self, num_frames, dt):
        self.initialize_animation()

        # Create the animation
        animation = FuncAnimation(self.fig, self.update_animation, frames=num_frames, interval=50, repeat=False)

        # Save the animation as an MP4 file
        animation.save('navierstoke_test.mp4', writer='ffmpeg', fps=30)
        
domain_size = 50
viscosity = 0.01
density = 1.0
dt = 0.01
solver = NavierStokesSolver(domain_size, viscosity, density)

# Animate the simulation for 100 frames with a time step of 0.01
solver.animate(num_frames=100, dt=dt)



#--------IN ONE DIMENSION
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class NavierStokesSolver1D:
    def __init__(self, domain_size, viscosity, density):
        self.domain_size = domain_size
        self.viscosity = viscosity
        self.density = density

        # Define grid
        self.grid_size = domain_size
        self.velocity_field = np.zeros(self.grid_size)
        self.pressure_field = np.zeros(self.grid_size)

        # Initialize with some initial conditions
        self.velocity_field[:] = 1.0  # Initial velocity

    def solve_step(self, dt):
        # Calculate velocity gradient
        du_dx = np.gradient(self.velocity_field, axis=0)

        # Calculate Laplacian of velocity field
        d2u_dx2 = np.gradient(du_dx, axis=0)

        # Calculate pressure change
        pressure_change = np.zeros_like(self.pressure_field)
        pressure_change[1:-1] = -self.density * (
            du_dx[1:-1] * (du_dx[1:-1] + d2u_dx2[1:-1])
        )

        # Update velocity field based on pressure change
        self.velocity_field[1:-1] += dt * (1 / self.density) * du_dx[1:-1] + dt * self.viscosity * d2u_dx2[1:-1]

        # Update pressure field
        self.pressure_field[1:-1] += dt * pressure_change[1:-1]

        # Enforce boundary conditions
        self.velocity_field[0] = 0  # Velocity is zero at the left boundary
        self.velocity_field[-1] = 0  # Velocity is zero at the right boundary

        self.pressure_field[0] = self.pressure_field[1]  # Pressure is constant at the left boundary
        self.pressure_field[-1] = self.pressure_field[-2]  # Pressure is constant at the right boundary

    def initialize_animation(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 4))

        # Initialize quiver plot for velocity field
        self.quiver = self.ax.quiver([], [], [], [], scale=20)
        self.ax.set_title('Velocity Field')
        self.ax.set_xlabel('X-axis')
        self.ax.set_ylabel('Velocity')

    def update_animation(self, frame):
        self.solve_step(dt)

        # Update quiver plot for velocity field
        positions = np.arange(0, self.domain_size)
        self.quiver.set_offsets(np.column_stack((positions[1:-1], np.zeros_like(self.velocity_field[1:-1]))))
        self.quiver.set_UVC(self.velocity_field[1:-1], np.zeros_like(self.velocity_field[1:-1]))

    def animate(self, num_frames, dt):
        self.initialize_animation()

        # Create the animation
        animation = FuncAnimation(self.fig, self.update_animation, frames=num_frames, interval=50, repeat=False)

        # Save the animation as an MP4 file
        animation.save('navierstoke_test.mp4', writer='ffmpeg', fps=30)

# Example Usage:
domain_size = 50
viscosity = 0.01
density = 1.0
dt = 0.01
num_frames = 100
solver = NavierStokesSolver1D(domain_size, viscosity, density)

# Animate the simulation for 100 frames with a time step of 0.01 and save as "navierstoke_test.mp4"
solver.animate(num_frames=num_frames, dt=dt)







import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class NavierStokesSolver1D:
    def __init__(self, domain_size, viscosity, density):
        self.domain_size = domain_size
        self.viscosity = viscosity
        self.density = density

        # Define grid
        self.grid_size = domain_size
        self.velocity_field = np.zeros(self.grid_size)
        self.pressure_field = np.zeros(self.grid_size)

        # Initialize with some initial conditions
        self.velocity_field[:] = 1.0  # Initial velocity

    def solve_step(self, dt):
        # Calculate velocity gradient
        du_dx = np.gradient(self.velocity_field, axis=0)

        # Calculate Laplacian of velocity field
        d2u_dx2 = np.gradient(du_dx, axis=0)

        # Calculate pressure change
        pressure_change = np.zeros_like(self.pressure_field)
        pressure_change[1:-1] = -self.density * (
            du_dx[1:-1] * (du_dx[1:-1] + d2u_dx2[1:-1])
        )

        # Update velocity field based on pressure change
        self.velocity_field[1:-1] += dt * (1 / self.density) * du_dx[1:-1] + dt * self.viscosity * d2u_dx2[1:-1]

        # Update pressure field
        self.pressure_field[1:-1] += dt * pressure_change[1:-1]

        # Enforce boundary conditions
        self.velocity_field[0] = 0  # Velocity is zero at the left boundary
        self.velocity_field[-1] = 0  # Velocity is zero at the right boundary

        self.pressure_field[0] = self.pressure_field[1]  # Pressure is constant at the left boundary
        self.pressure_field[-1] = self.pressure_field[-2]  # Pressure is constant at the right boundary

    def initialize_animation(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 4))

        # Initialize quiver plot for velocity field
        self.quiver = self.ax.quiver([], [], [], [], scale=20)
        self.ax.set_title('Velocity Field')
        self.ax.set_xlabel('X-axis')
        self.ax.set_ylabel('Velocity')

    def update_animation(self, frame):
        self.solve_step(dt)

        # Update quiver plot for velocity field
        positions = np.arange(0, self.domain_size)
        self.quiver.set_offsets(np.column_stack((positions[1:-1], np.zeros_like(self.velocity_field[1:-1]))))
        self.quiver.set_UVC(self.velocity_field[1:-1], np.zeros_like(self.velocity_field[1:-1]), np.zeros_like(self.velocity_field[1:-1]))

    def animate(self, num_frames, dt):
        self.initialize_animation()

        # Create the animation
        animation = FuncAnimation(self.fig, self.update_animation, frames=num_frames, interval=50, repeat=False)

        # Save the animation as an MP4 file
        animation.save('navierstoke_test.mp4', writer='ffmpeg', fps=30)

# Example Usage:
domain_size = 50
viscosity = 0.01
density = 1.0
dt = 0.01
num_frames = 100
solver = NavierStokesSolver1D(domain_size, viscosity, density)

# Animate the simulation for 100 frames with a time step of 0.01 and save as "navierstoke_test.mp4"
solver.animate(num_frames=num_frames, dt=dt)



#-------woirkable------
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

class HeatTransferSimulation:
    def __init__(self, initial_temperature, heat_transfer_coefficient, time_step, num_dimensions):
        self.temperature = np.full(num_dimensions, initial_temperature)
        self.heat_transfer_coefficient = heat_transfer_coefficient
        self.time_step = time_step
        self.num_dimensions = num_dimensions
        self.temperature_history = [self.temperature.copy()]

    def update_temperature(self, external_temperature):
        temperature_difference = external_temperature - self.temperature
        heat_transfer = self.heat_transfer_coefficient * temperature_difference * self.time_step
        self.temperature += heat_transfer

    def run_simulation(self, external_temperature, num_steps):
        for _ in range(num_steps):
            self.update_temperature(external_temperature)
            self.temperature_history.append(self.temperature.copy())

    def animate_temperature_over_time(self):
        fig, ax = plt.subplots()
        dimensions = range(1, self.num_dimensions + 1)
        line, = ax.plot(dimensions, self.temperature_history[0])

        def update(frame):
            line.set_ydata(self.temperature_history[frame])
            ax.set_title(f'Temperature Evolution - Time Step {frame + 1}')
            return line,

        def init():
            ax.set_ylim(0, 150)  # Set y-axis limits to unlimited initially
            return line,

        anim = FuncAnimation(fig, update, frames=len(self.temperature_history), init_func=init, interval=500, blit=True)
        plt.xlabel("Dimension")
        plt.ylabel("Temperature")

        # Save animation as a GIF
        anim.save("heated_simm.gif", writer="pillow", fps=1)

# Example usage
initial_temperature = 100.0
heat_transfer_coefficient = 0.02
time_step = 1
num_dimensions = 5
num_steps = 60

simulation = HeatTransferSimulation(initial_temperature, heat_transfer_coefficient, time_step, num_dimensions)
simulation.run_simulation(external_temperature=25.0, num_steps=num_steps)
simulation.animate_temperature_over_time()





#------real one------------


import numpy as np

import matplotlib.pyplot as plt

#we have node as continous and each node is we meausre the temperature

#length and alpha

a=110 #copper
length = 50 #mm
time= 4 #what time we look at that second
nodes=10 #this is delta x

#WE HAVE ALSO DELTA T for that

#**smaller delta x better precision

#delta t can quickly

dx=length/nodes
#dx shode5
#dt shode
dt=0.5*dx**2/a
#yani har 0.11 sanie


t_nodes=int(time/dt)
#t_nodes shode 35
#yani 35 ta baze zamani darim

#boundary
#temperature of boundary at any time is certain temperature:
    #but other is initial at t=0

u=np.zeros(nodes)+20
#10 ta 20 misaze ke hame 20 hastan
#20 darajan hame
u[0]=100
u[-1]=100
#ina hamishe 100 astan

#simulation

counter=0
#t azamani k zire zaman bashim
while counter<time:
    #we need seperate copy of distribution
    w=u.copy()
    
    #baraye node 1 ta yeki monde be akahr
    #fromole
    #Ti+1,i = Ti,i + alpha dar 
    #yani t jadid ba gahblie ine formolesh
    
    for i in range(1,nodes-1):
        u[i]=dt*a*(w[i-1]-2*w[i]+w[i+1])/dx**2+w[i]
        
    counter+=dt
    
    #print('t:',counter,' [s]   ', 'aVERAGE tEMP: ',np.average(u))
    print(" t: {:.3f} [s], Average Temperature: {:.2f} C".format(counter,np.average(u)))



#======================================================
#visual
#copper

import numpy as np

import matplotlib.pyplot as plt
#we have node as continous and each node is we meausre the temperature
#length and alpha
a=0.34 #copper

#copper 110 
#alminium 97
length = 50 #mm
time= 120 #what time we look at that second
nodes=10 #this is delta x
#WE HAVE ALSO DELTA T for that
#**smaller delta x better precision
#delta t can quickly
dx=length/nodes
#dx shode5
#dt shode
#dt=0.5*dx**2/a
dt=0.2
#yani har 0.11 sanie
t_nodes=int(time/dt)
#t_nodes shode 35
#yani 35 ta baze zamani darim
#boundary
#temperature of boundary at any time is certain temperature:
    #but other is initial at t=0
u=np.zeros(nodes)+20
#10 ta 20 misaze ke hame 20 hastan
#20 darajan hame
u[0]=100
u[-1]=100
#ina hamishe 100 astan
fig, axis=plt.subplots()
pcm = axis.pcolormesh([u],cmap=plt.cm.jet,vmin=0,vmax=130)
plt.colorbar(pcm,ax=axis)
axis.set_ylim(-2,3) #optional
#simulation
counter=0
#t azamani k zire zaman bashim
while counter<time:
    #we need seperate copy of distribution
    w=u.copy()
    
    #baraye node 1 ta yeki monde be akahr
    #fromole
    #Ti+1,i = Ti,i + alpha dar 
    #yani t jadid ba gahblie ine formolesh
    
    for i in range(1,nodes-1):
        u[i]=dt*a*(w[i-1]-2*w[i]+w[i+1])/dx**2+w[i]
        
    counter+=dt
    
    #print('t:',counter,' [s]   ', 'aVERAGE tEMP: ',np.average(u))
    #print(" t: {:.3f} [s], Average Temperature: {:.2f} C".format(counter,np.average(u)))
    
    #update
    pcm.set_array([u])
    axis.set_title('distribution at t: {:.3f}'.format(counter))
    plt.pause(0.0001)




#glas--------------------------------

import numpy as np

import matplotlib.pyplot as plt
#we have node as continous and each node is we meausre the temperature
#length and alpha
a=0.34 #copper

#copper 110 
#alminium 97
length = 50 #mm
time= 120 #what time we look at that second
nodes=10 #this is delta x
#WE HAVE ALSO DELTA T for that
#**smaller delta x better precision
#delta t can quickly
dx=length/nodes
#dx shode5
#dt shode
#dt=0.5*dx**2/a
dt=0.2
#yani har 0.11 sanie
t_nodes=int(time/dt)
#t_nodes shode 35
#yani 35 ta baze zamani darim
#boundary
#temperature of boundary at any time is certain temperature:
    #but other is initial at t=0
u=np.zeros(nodes)+20
#10 ta 20 misaze ke hame 20 hastan
#20 darajan hame
u[0]=100
u[-1]=100
#ina hamishe 100 astan
fig, axis=plt.subplots()
pcm = axis.pcolormesh([u],cmap=plt.cm.jet,vmin=0,vmax=130)
plt.colorbar(pcm,ax=axis)
axis.set_ylim(-2,3) #optional
#simulation
counter=0
#t azamani k zire zaman bashim
while counter<time:
    #we need seperate copy of distribution
    w=u.copy()
    
    #baraye node 1 ta yeki monde be akahr
    #fromole
    #Ti+1,i = Ti,i + alpha dar 
    #yani t jadid ba gahblie ine formolesh
    
    for i in range(1,nodes-1):
        u[i]=dt*a*(w[i-1]-2*w[i]+w[i+1])/dx**2+w[i]
        
    counter+=dt
    
    #print('t:',counter,' [s]   ', 'aVERAGE tEMP: ',np.average(u))
    #print(" t: {:.3f} [s], Average Temperature: {:.2f} C".format(counter,np.average(u)))
    
    #update
    pcm.set_array([u])
    axis.set_title('distribution at t: {:.3f}'.format(counter))
    plt.pause(0.0001)







def d1_heat_conductivity(material,left_temp,right_temp,initial_temp,leng,timee,step='def'):
    
    if material=='silver':
        a=149
    if material=='gold':
        a=127
    if material=='copper':
        a=113
    if material=='alminium':
        a=97.5
    if material=='iron':
        a=22.8
    if material=='mercury':
        a=4.7
    if material=='marble':
        a=1.2
    if material=='ice':
        a=1.2
    if material=='concrtee':
        a=0.75
    if material=='brick':
        a=0.52
    if material=='glass':
        a=0.34
    if material=='wood':
        a=0.13

    length = leng #mm
    time= timee #what time we look at that second
    nodes=10 #this is delta x
    #WE HAVE ALSO DELTA T for that
    #**smaller delta x better precision
    #delta t can quickly
    dx=length/nodes
    #dx shode5
    #dt shode
    if a>80:
        dt=0.5*dx**2/a
    else:
        dt=0.2
    if step=='def':
        dtt=0.0001
    else:
        dtt=dt

    #yani har 0.11 sanie
    t_nodes=int(time/dt)
    #t_nodes shode 35
    #yani 35 ta baze zamani darim
    #boundary
    #temperature of boundary at any time is certain temperature:
        #but other is initial at t=0
    u=np.zeros(nodes)+ initial_temp
    #10 ta 20 misaze ke hame 20 hastan
    #20 darajan hame
    u[0]=left_temp
    u[-1]=right_temp
    #ina hamishe 100 astan
    fig, axis=plt.subplots()
    pcm = axis.pcolormesh([u],cmap=plt.cm.jet,vmin=0,vmax=130)
    plt.colorbar(pcm,ax=axis)
    axis.set_ylim(-2,3) #optional
    #simulation
    counter=0
    #t azamani k zire zaman bashim
    while counter<time:
        #we need seperate copy of distribution
        w=u.copy()
        
        #baraye node 1 ta yeki monde be akahr
        #fromole
        #Ti+1,i = Ti,i + alpha dar 
        #yani t jadid ba gahblie ine formolesh
        
        for i in range(1,nodes-1):
            u[i]=dt*a*(w[i-1]-2*w[i]+w[i+1])/dx**2+w[i]
            
        counter+=dt
        
        #print('t:',counter,' [s]   ', 'aVERAGE tEMP: ',np.average(u))
        #print(" t: {:.3f} [s], Average Temperature: {:.2f} C".format(counter,np.average(u)))
        
        #update
        pcm.set_array([u])
        axis.set_title('distribution at t: {:.3f}'.format(counter))
        plt.pause(dtt)
    







d1_heat_conductivity('gold', 150, 90, 25, 0.1 , 40,step='def')


import numpy as np
import matplotlib.pyplot as plt

def d2_heat_conductivity(material,upper_temp,lower_temp,left_temp,right_temp,initial_temp,leng,timee,nod,step='def'):

    if material=='silver':
        a=149
    if material=='gold':
        a=127
    if material=='copper':
        a=113
    if material=='alminium':
        a=97.5
    if material=='iron':
        a=22.8
    if material=='mercury':
        a=4.7
    if material=='marble':
        a=1.2
    if material=='ice':
        a=1.2
    if material=='concrtee':
        a=0.75
    if material=='brick':
        a=0.52
    if material=='glass':
        a=0.34
    if material=='wood':
        a=0.13

    length = leng #mm
    time= timee #what time we look at that second
    nodes=nod #this is delta x
    #WE HAVE ALSO DELTA T for that
    #**smaller delta x better precision
    #delta t can quickly
    dx=length/nodes
    dy=length/nodes
    
    #dx shode5
    #dt shode
    if a>80:
        dt=min(dx**2/(4*a),dy**2/(4*a))
    else:
        dt=0.2
        
    if step=='def':
        #dtt=dt
        dtt=10**-12
    else:
        dtt=dt

    #yani har 0.11 sanie
    t_nodes=int(time/dt)
    #t_nodes shode 35
    #yani 35 ta baze zamani darim
    #boundary
    #temperature of boundary at any time is certain temperature:
        #but other is initial at t=0
    #**ta alan ye dy ezafe krdim
    #hala 2d mikonim
    u=np.zeros((nodes,nodes))+ initial_temp
    #10 ta 20 misaze ke hame 20 hastan
    #20 darajan hame
    
    #mishe chartaee kard
    u[0,:]=lower_temp
    u[-1,:]=upper_temp
    u[:,0]=left_temp
    u[:,-1]=right_temp
    
    
    
    #ina hamishe 100 astan
    fig, axis=plt.subplots()
    #remove [] from u
    pcm = axis.pcolormesh(u,cmap=plt.cm.jet,vmin=0,vmax=130)
    plt.colorbar(pcm,ax=axis)
    #remove y_lim
    
    
    #simulation
    counter=0
    #t azamani k zire zaman bashim
    while counter<time:
        #we need seperate copy of distribution
        w=u.copy()
        
        #baraye node 1 ta yeki monde be akahr
        #fromole
        #Ti+1,i = Ti,i + alpha dar 
        #yani t jadid ba gahblie ine formolesh
        
        for i in range(1,nodes-1):
            #another j
            for j in range(1,nodes-1):
                
                #u[i]=dt*a*(w[i-1]-2*w[i]+w[i+1])/dx**2+w[i]
                
                dd_ux=(w[i-1,j]-2*w[i,j]+w[i+1,j])/dx**2
                dd_uy=(w[i,j-1]-2*w[i,j]+w[i,j+1])/dy**2
                
                u[i,j]=dt*a*(dd_ux+dd_uy)+w[i,j]
            
        counter+=dt
        
        #print('t:',counter,' [s]   ', 'aVERAGE tEMP: ',np.average(u))
        #print(" t: {:.3f} [s], Average Temperature: {:.2f} C".format(counter,np.average(u)))
        
        #update
        #remove the bracket
        pcm.set_array(u)
        axis.set_title('distribution at t: {:.3f}'.format(counter))
        plt.pause(dtt)
    


d2_heat_conductivity('gold', 150,25,25,150, 25, 50 , 60,40,step='def')

d2_heat_conductivity('alminium', 150,25,25,150, 25, 50 , 60,40,step='def')


d2_heat_conductivity('glass', 150,25,25,150, 25, 50 , 60,40,step='def')


d2_heat_conductivity('glass', 100,100,25,25, 25, 50 , 30,40,step='deffg')


#then we can creat a class to do that

#-------NAVIE STOKE
#linear
import numpy 
from matplotlib import pyplot as plt

nx=41 #411 to 81
dx=2/(nx-1)
nt=25 #time
dt=0.025
c=1  #wqave speed of c=1

u=numpy.ones(nx)
u[int(.5/dx):int(1/dx+1)]=2 
print(u)
plt.plot(numpy.linspace(0, 2,nx),u)

un=numpy.ones(nx) #initiate temporary

for n in range(nt):
    un=u.copy()
    for i in range(1,nx):
        u[i]=un[i]-c*dt/dx*(un[i]-un[i-1])
    plt.plot(numpy.linspace(0,2,nx),u);
    plt.pause(0.1)
print(u)
plt.show()


#non-linear
import numpy 
from matplotlib import pyplot as plt

nx=41 #41 to 81
dx=2/(nx-1)
nt=25 #time
dt=0.025
c=1  #wqave speed of c=1

u=numpy.ones(nx)
u[int(.5/dx):int(1/dx+1)]=2 
print(u)
plt.plot(numpy.linspace(0, 2,nx),u)

un=numpy.ones(nx) #initiate temporary

for n in range(nt):
    un=u.copy()
    for i in range(1,nx):
        #just this
        u[i]=un[i]-un[i]*dt/dx*(un[i]-un[i-1])
    plt.plot(numpy.linspace(0,2,nx),u);
    plt.pause(0.01)
print(u)
plt.show()



#https://www.youtube.com/watch?v=9XHkaHkRAOg&list=PLE4jpqcRJiBpODH_ksfgJmKsSc9-CN3_A&index=3
#if you chnage the nx it change
#difgfusion

import numpy              			   
from matplotlib import pyplot as plt 

# Initial Conditions
nx = 41
dx = 2 / (nx - 1)
nt = 20    #the number of timesteps we want to calculate
nu = 0.3   #the value of viscosity
sigma = .2 #sigma is a parameter, we'll learn more about it later
dt = sigma * dx**2 /nu #dt is defined using sigma ... more later!
print(dt)

u = numpy.ones(nx)      #a numpy array with nx elements all equal to 1.
u[int(.5 / dx):int(1 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s

# Calculation
un = numpy.ones(nx) #our placeholder array, un, to advance the solution in time
for n in range(nt):  #iterate through time
    un = u.copy() ##copy the existing values of u into un
    for i in range(1, nx - 1):
        u[i] = un[i] + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])
        plt.pause(0.01)
        
    plt.plot(numpy.linspace(0, 2, nx), u);
plt.show();


#2d
from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots
import numpy
from matplotlib import pyplot
from matplotlib import cm 
import time

###variable declarations
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

#sakhtane x o y
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx)) ##create a 1xn vector of 1's
un = numpy.ones((ny, nx)) ##

###Assign initial conditions
##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2 

###Plot Initial Condition
##the figsize parameter can be used to produce different sized images
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')     
X, Y = numpy.meshgrid(x, y)                            
surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
pyplot.show();


# Calc with Loops
t = time.time()
u = numpy.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

for n in range(nt + 1): ##loop across number of time steps
    un = u.copy()
    row, col = u.shape
    for j in range(1, row):
        for i in range(1, col):
            u[j, i] = (un[j, i] - (c * dt / dx * (un[j, i] - un[j, i - 1])) -
                                  (c * dt / dy * (un[j, i] - un[j - 1, i])))
            u[0, :] = 1
            u[-1, :] = 1
            u[:, 0] = 1
            u[:, -1] = 1

fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
pyplot.show();
print('The Calculation with loops took: ' + str('{0:.2f}'.format(time.time() - t)) + 's')




# Calc with Arrays
t = time.time()
u = numpy.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

for n in range(nt + 1): ##loop across number of time steps
    un = u.copy()
    u[1:, 1:] = (un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, :-1])) -
                              (c * dt / dy * (un[1:, 1:] - un[:-1, 1:])))
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
pyplot.show();
print('The Calculation with arrays took: ' + str('{0:.2f}'.format(time.time() - t)) + 's')




#======
import numpy
from matplotlib import pyplot
from matplotlib import cm 
from mpl_toolkits.mplot3d import Axes3D
import time
### variable declarations
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = 0.2
dt = sigma * dx

# sakhtane x o y
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))  ## create a 1xn vector of 1's
un = numpy.ones((ny, nx))

### Assign initial conditions
## set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(0.5 / dy):int(1 / dy + 1), int(0.5 / dx):int(1 / dx + 1)] = 2 

### Plot Initial Condition
## the figsize parameter can be used to produce different sized images
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.add_subplot(111, projection='3d')  # Use add_subplot instead of gca
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
pyplot.show()

'''



#**
'''
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')     
X, Y = numpy.meshgrid(x, y)                            
surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
pyplot.show();
'''




'''
# Calc with Loops
t = time.time()
u = numpy.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

for n in range(nt + 1): ##loop across number of time steps
    un = u.copy()
    row, col = u.shape
    for j in range(1, row):
        for i in range(1, col):
            u[j, i] = (un[j, i] - (c * dt / dx * (un[j, i] - un[j, i - 1])) -
                                  (c * dt / dy * (un[j, i] - un[j - 1, i])))
            u[0, :] = 1
            u[-1, :] = 1
            u[:, 0] = 1
            u[:, -1] = 1
            
            fig = pyplot.figure(figsize=(11, 7), dpi=100)
            #ax = fig.gca(projection='3d')
            ax = fig.add_subplot(111, projection='3d') 
            surf2 = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
            pyplot.show();

            
pyplot.show()


fig = pyplot.figure(figsize=(11, 7), dpi=100)
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(111, projection='3d') 
surf2 = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
pyplot.show();
print('The Calculation with loops took: ' + str('{0:.2f}'.format(time.time() - t)) + 's')








#maybe animation-----------------
import numpy
from matplotlib import pyplot
from matplotlib import cm 
from mpl_toolkits.mplot3d import Axes3D
import time

### variable declarations
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = 0.2
dt = sigma * dx

# sakhtane x o y
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))  ## create a 1xn vector of 1's
un = numpy.ones((ny, nx))

### Assign initial conditions
## set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(0.5 / dy):int(1 / dy + 1), int(0.5 / dx):int(1 / dx + 1)] = 2 

### Plot Initial Condition
## the figsize parameter can be used to produce different sized images

fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.add_subplot(111, projection='3d')  # Use add_subplot instead of gca
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
pyplot.show()







# Create the figure and subplot outside the loop
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.add_subplot(111, projection='3d') 

# Calc with Loops
t = time.time()
u = numpy.ones((ny, nx))
u[int(0.5 / dy):int(1 / dy + 1), int(0.5 / dx):int(1 / dx + 1)] = 2

for n in range(nt + 1): ##loop across number of time steps
    un = u.copy()
    row, col = u.shape
    for j in range(1, row):
        for i in range(1, col):
            u[j, i] = (un[j, i] - (c * dt / dx * (un[j, i] - un[j, i - 1])) -
                                  (c * dt / dy * (un[j, i] - un[j - 1, i])))
            u[0, :] = 1
            u[-1, :] = 1
            u[:, 0] = 1
            u[:, -1] = 1
            
    ax.clear()  # Clear the previous plot
    surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
    pyplot.pause(0.1)  # Pause for a short time to create animation effect

pyplot.show()



#molecular dynamic
#https://github.com/PolymerTheory/MDFromScratch/blob/main/MDFS.ipynb


import numpy as np

n=5 # 5 molecules
D=3 # 3 dimention
LL=20.0 #in Angstroms
dt=0.0001 #In ps
dt=0.07

#arrays of variables
r=np.random.rand(n,D)*LL #initialize all random positions
v=np.random.rand(n,D)-0.5 #initialize with random velocity
'''


#**
'''
def update(r,v,dt):
    newr=r+dt*v
    if BC==0:
        newr=newr%L
        newv=1.0*v
    if BC==1:
        newr,newv=reflectBC(newr,v)
    return newr,newv

'''


'''
def update(r,v,dt):
    newr=r+dt*v
    return newr


for i in range(100):
    r=update(r,v,dt)
    print(i)
    print(r)

from mpl_toolkits import mplot3d
import pylab as pl
from IPython import display
import matplotlib.pyplot as plt

#dobare v or
L=np.zeros([D])+LL

r=LL*np.random.rand(n,D)
v=100.0*(np.random.rand(n,D)-0.5)


fig = plt.figure()
for i in range(100):
    r=update(r,v,dt)
    #ax=plt.axes(projection='3d')
    ax = fig.add_subplot(111, projection='3d') 
    ax.set_xlim3d(0,L[0])
    ax.set_ylim3d(0,L[1])
    ax.set_zlim3d(0,L[2])
    ax.scatter3D([item[0] for item in r],[item[1] for item in r],[item[2] for item in r])
    #display.clear_output(wait=True)
   # display.display(pl.gcf())
    plt.pause(dt)
plt.show()
#bayad 3d ish kard


#vaghena dt bayad yekari krd vasash
BC=1 # 0 for periodic, 1 for reflecting

#for boundary
def update(r,v,dt):
    newr=r+dt*v
    if BC==0:
        newr=newr%L
        newv=1.0*v
    return newr


#make particles reflect off boundaries
def reflectBC(r,v):
    newv = 1.0*v
    newr = 1.0*r
    for i in range(n):
        for j in range(D):
            #for lowere
            if newr[i][j]<0:
                newr[i][j]= -newr[i][j]
                newv[i][j]=abs(v[i][j])
            #for upper
            if newr[i][j]>L[j]:
                newr[i][j]= 2.0*L[j]-newr[i][j]
                newv[i][j]=-abs(v[i][j])
    return newr,newv

def update(r,v,dt):
    newr=r+dt*v
    if BC==0:
        newr=newr%L
        newv=1.0*v
    if BC==1:
        newr,newv=reflectBC(newr,v)
    return newr,newv



#visualizesh with avodi

def dump(r,t):
    fname=outdir+"/t"+str(t)+".dump"
    f=open(fname,"w")
    f.write("ITEM: TIMESTEP\n")
    f.write(str(t)+"\n") #time step
    f.write("ITEM: NUMBER OF ATOMS\n")
    f.write(str(len(r))+"\n") # number of atoms
    #**maybe change
    f.write("ITEM: BOX BOUNDS pp pp pp\n") #pp = periodic BCs
    f.write("0 "+str(L[0])+"\n")
    f.write("0 "+str(L[1])+"\n")
    f.write("0 "+str(L[2])+"\n")
    f.write("ITEM: ATOMS id mol type x y z\n")
    for i in range(len(r)):
        f.write(str(i)+" "+str(mols[i])+" "+str(tp[i])+" "+str(r[i][0])+" "+str(r[i][1])+" "+str(r[i][2])+"\n")
    f.close
    

    
for i in range(1000):
    r,v=update(r,v,dt)
    dump(r,int(i))

#tehn add 1000 line to one code
#then go to avodio nd then the file.dump
#for each dumb file


#you can play all .dump for video


#so far non-contacting
#if we have we have Rij -->leonard jones
#rij=[xi-xj + yi-yj + ... radical--> this is distance between two atom

#P is pressure
#Pij=..... that has rij and sigma on it for two particle

#potantiel energy (Ui) is sum of all P that one i has with other j (not itself) zigmaP
#so the Fi--> - gradient Ui


#so dP/dX=-24...
#FOR Y AND Z IS THE SAME
#WE can have the matrixes --> for gradiant Pij=-24 rij sima * xij yij zij



#leonard jones potentisal
eps=1.0
sig=1.0

#particle mass=1.0
mass=1.0

#this is calculate the P = 4 sigma sigmma/r ...
#tp=


def LJpot(r,i,sig,eps):
    drv=r-r[i]
    drv=np.delete(drv,i,0)
    dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
    r6=(sig/np.array(dr))**6
    r12=(sig/np.array(dr))**12
    LJP=4.0*eps*sum(r12-r6)
    return LJP
    
'''




#**
'''
def LJpot(r,i,sigg,epss):
    sg=np.delete(np.array([sigg[tp[j]] for j in range(n)]),i)
    ep=np.delete(np.array([epss[tp[j]] for j in range(n)]),i)
    for ii in range(n): #ignore atoms in the same molecule
        if mols[i]==mols[ii]:
            ep[ii]=0
    drv=r-r[i] #distance in each dimension
    drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
    dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
    r6 =(sg/np.array(dr))**6
    r12=(sg/np.array(dr))**12
    LJP=4.0*eps*sum(ep*r6-ep*r12)
    return LJP
'''




'''
def dLJp(r,u,sig,eps):
    drv=r-r[i] #distance in each dimension
    drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
    dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
    r8 =    (sig**6)*(1.0/dr)**8
    r14=2.0*(sig**12)*(1.0/dr)**14
    r814=r14-r8
    r814v=np.transpose(np.transpose(drv)*r814)
    r814vs=np.sum(r814v,axis=0)
    dLJP=24.0*eps*(r814vs)
    return dLJP




#this is calculate the gradiat --> delta P=-24... xij yij zij
#Gradient of Lennard-Jones potential
'''


#**
'''
def dLJp(r,i,sigl,epsl,bnds):
    sg=np.delete(np.array([sigl[tp[j]] for j in range(n)]),i)
    ep=np.array([epsl[tp[j]] for j in range(n)])
    for ii in range(n): #ignore atoms in the same molecule
        if mols[i]==mols[ii]:
            ep[ii]=0
    ep=np.delete(ep,i)
    drv=r-r[i] #distance in each dimension
    drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
    dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
    r8 =    ep*(sg**6)*(1.0/np.array(dr))**8
    r14=2.0*ep*(sg**12)*(1.0/np.array(dr))**14
    r8v =np.transpose(np.transpose(drv)*r8)
    r14v=np.transpose(np.transpose(drv)*r14)
    r8vs =np.sum(r8v,axis=0)
    r14vs=np.sum(r14v,axis=0)
    dLJP=24.0*(r14vs-r8vs)
    return dLJP

'''



'''

#both def getting the i

#for each I we can getthe P and all Gradiant P

#Update subroutine
#Update subroutines
def updatev(r,v,dt,sig,eps):
    #calculate acceleration:
    F=-np.array([dLJp(r,i,sig,eps) for i in range(n)]) #LJ
    a=F/mass
    newv=v+a*dt
    return newv,a

def update(r,v,dt):
    newr=r+dt*v
    if BC==0:
        newr=newr%L
        newv=1.0*v
    if BC==1:
        newr,newv=reflectBC(newr,v)
    return newr,newv

#har noghte choon y shetabi dare az baghie migire taghir mikone na random dg
#pas b hesabe v ham r ham taghir mikone


'''
'''




#**
def updatev(r,v,dt,sigg,epss):
    #calculate acceleration:
    F=-np.array([dLJp(r,i,sigg[tp[i]],epss[tp[i]],bnd) for i in range(n)]) #LJ
    F=F-dBEpot(r,bnd) #Bonds
    F=F-dBA(r,angs) #Bonds angles
    F=F-np.array([coul(r,i,chrg) for i in range(n)]) #Coulomb
    a=np.transpose(np.transpose(F)/mm) #Force->acceleration
    #update velocity
    newv=v+dt*a
    return newv,a
'''

#tp for periodic

#=========================================
#=========================================
#=========================================
#=========================================
#========example=======
#=========================================
#=========================================
#=========================================

'''
#you can start with two particcle

from mpl_toolkits import mplot3d
import pylab as pl
from IPython import display
import matplotlib.pyplot as plt
import numpy as np

n=2 # 5 molecules
D=3 # 3 dimention
LL=10.0 #in Angstroms
#dt=0.07 #in Ps
mass=1

# dota particle hastan ke 3 ta fasele beyneshone
#hey beham mikhoran barmirgdn
#arrays of variables
#r=np.random.rand(n,D)*LL r ro tarif mikonim
r=np.zeros((n,D))  #for avoiding the error --> NameError: name 'r' is not defined
r[0][0]=LL/2 + 1.5
r[0][1]=LL/2 
r[0][2]=LL/2 


r[1][0]=LL/2 - 1.5
r[1][1]=LL/2 
r[1][2]=LL/2 

v=np.zeros((n,D))  #initialize with 0


#dobare v or
L=np.zeros([D])+LL #in chie???
#r=LL*np.random.rand(n,D) -->niazi nis darim
#v=100.0*(np.random.rand(n,D)-0.5) niazi nis darim
sigg=1.0
epss=1.0
dtt=0.07
dt=0.07
BC=1

#dar v mikhaymesh
def dLJp(r,u,sig,eps):
    drv=r-r[i] #distance in each dimension
    drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
    dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
    r8 =(sig**6)*(1.0/np.array(dr))**8
    r14=2.0*(sig**12)*(1.0/np.array(dr))**14
    r814=r14-r8
    r814v=np.transpose(np.transpose(drv)*r814)
    r814vs=np.sum(r814v,axis=0)
    dLJP=24.0*eps*(r814vs)
    return dLJP

sig=1.0
eps=1.0
i=0



def updatev(r,v,dt,sig,eps):
    #calculate acceleration:
    F=-np.array([dLJp(r,i,sig,eps) for i in range(n)]) #LJ
    a=F/mass
    newv=v+a*dt
    return newv,a
    #return newv ya in ya brim paeen


def reflectBC(r,v):
    newv = 1.0*v
    newr = 1.0*r
    for i in range(n):
        for j in range(D):
            #for lowere
            if newr[i][j]<0:
                newr[i][j]= -newr[i][j]
                newv[i][j]=abs(v[i][j])
            #for upper
            if newr[i][j]>L[j]:
                newr[i][j]= 2.0*L[j]-newr[i][j]
                newv[i][j]=-abs(v[i][j])
    return newr,newv


def update(r,v,dt):
    newr=r+dt*v
    if BC==0:
        newr=newr%L
        newv=1.0*v
    if BC==1:
        newr,newv=reflectBC(newr,v)
    return newr,newv



i=0

fig = plt.figure()
for i in range(100):
    #inja bejaye fght update e r , v ham update mikonim
    #ama aval v ro update kon
    #vn yani new
    v=updatev(r,v,dtt,sigg,epss)[0]
    r=update(r,v,dt)[0]
    #????? v e chijori jadido mziaare o ina
    
    #ax=plt.axes(projection='3d')
    ax = fig.add_subplot(111, projection='3d') 
    ax.set_xlim3d(0,L[0])
    ax.set_ylim3d(0,L[1])
    ax.set_zlim3d(0,L[2])
    ax.scatter3D([item[0] for item in r],[item[1] for item in r],[item[2] for item in r])
    #display.clear_output(wait=True)
   # display.display(pl.gcf())
    plt.pause(dt)
plt.show()

#hamin nshon mide ba yekam moshkel

'''
