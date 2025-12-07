import numpy as np
import matplotlib.pyplot as plt

print("Kinetics Simulator Module Loaded")
class ReactionKineticsSimulator:
    def __init__(self, C0, k, t_end, dt, reaction_type="first_order", k_rev=None):
        self.C0 = C0
        self.k = k
        self.k_rev = k_rev
        self.t_end = t_end
        self.dt = dt
        self.reaction_type = reaction_type

        self.t = np.arange(0, t_end + dt, dt)
        self.C = np.zeros_like(self.t)
        self.C[0] = C0

    def rate(self, C):
        #Returns dC/dt depending on reaction type.
        if self.reaction_type == "first_order":
            return -self.k * C

        elif self.reaction_type == "second_order":
            return -self.k * C**2

        elif self.reaction_type == "reversible":
            return -self.k * C + self.k_rev * (self.C0 - C)

        else:
            raise ValueError("Invalid reaction type")

    def run(self):
        #Euler's Method Simulation
        for i in range(1, len(self.t)):
            dCdt = self.rate(self.C[i - 1])
            self.C[i] = self.C[i - 1] + dCdt * self.dt

        return self.t, self.C

    def plot(self):
        plt.plot(self.t, self.C, linewidth=2)
        plt.title(f"{self.reaction_type.replace('_', ' ').title()} Reaction")
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.grid(True)
        plt.show()

if __name__ == "__main__":
    print("Select reaction type:")
    print("1. First Order")
    print("2. Second Order")
    print("3. Reversible Reaction")

    choice = int(input("Enter choice (1/2/3): "))

    if choice == 1:
        r_type = "first_order"
        k_rev = None
    elif choice == 2:
        r_type = "second_order"
        k_rev = None
    elif choice == 3:
        r_type = "reversible"
        k_rev = float(input("Enter reverse rate constant k_reverse: "))
    else:
        raise ValueError("Invalid selection")

    C0 = float(input("Enter initial concentration (mol/L): "))
    k = float(input("Enter forward rate constant k: "))
    t_end = float(input("Enter total time of simulation: "))
    dt = float(input("Enter time step: "))

    sim = ReactionKineticsSimulator(
        C0=C0,
        k=k,
        k_rev=k_rev,
        t_end=t_end,
        dt=dt,
        reaction_type=r_type
    )

    t, C = sim.run()
    sim.plot()
