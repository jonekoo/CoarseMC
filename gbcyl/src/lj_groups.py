#-*- coding: utf-8 -*-
from abc import ABCMeta, abstractmethod
import copy
from random import random as rr


class ParticleGroup:
    __metaclass__ = ABCMeta

    def __init__(self, particles, groupID):
        self.particles = particles
        self.groupID = groupID
        self.interactions = []

    @abstractmethod
    def move(self):
        pass

    @abstractmethod
    def map(self, func):
        """Apply operation func to all particles."""
        pass

    @abstractmethod
    def addInteraction(ia):
        pass

    @abstractmethod
    def pairReduce(self, ia, cutoff=None):
        """Loop through all pairs of particles within cutoff."""
        pass

    @abstractmethod
    def reduce(self, ia, pos, cutoff=None):
        """Loop through all neighbours of pos (within cutoff)."""
        pass

class SphereGroup(ParticleGroup):
    
    def addInteraction(self, ia):
        self.interactions.append(ia)
    
    def move(self, rule):
        for p in self.particles:
            p0 = copy.deepcopy(p)
            e0 = self.energy(p)
            p.move()
            e1 = self.energy(p)
            if not rule.accept(e0, e1):
                p = p0

    def pairReduce(self, potential, cutoff=None):
        sum = 0.0
        for i, p in enumerate(self.particles):
            for o in self.particles[:i]:
                sum += potential(p, o)
        return sum

    def energy(self, particle):
        energy = 0.0
        for ia in self.interactions:
            energy = energy + ia.energyCallback(particle)
        return energy

    def reduce(self, potential, pos, cutoff):
        # We could leave out the pos and cutoff and include them in the
        # potential if it were an object. Potential could also contain
        # the other particle if needed. Cutoff checking could be left out
        # altogether and assume nearest neighbour reduce unless a cutoff
        # is explicitly given. Potential could also take a Sphere
        # as an argument which would simplify this function. It may not be
        # optimal for performance, though.
        sum = 0
        for p in self.particles:
            if ((p.x - pos[0])**2 + (p.y - pos[1])**2 + \
                    (p.z - pos[2])**2) < cutoff**2:
                # Jälkimmäinen argumentti callback-kutsussa ehkä turha
                sum = sum + potential.value(p.position())
        return sum

    # Let's now define a function which would allow us to modify all the
    # particles in the group. Here the iteration pattern is simple but more
    # complex patterns would be equally doable. Even move could basically
    # be defined as a map operation.
    def map(self, func):
        for p in self.particles:
            func(p)


class Sphere:

    max1D = 0.1

    def __init__(self, x=0.0, y=0.0, z=0.0, max1D = 0.1):
        self.x = x
        self.y = y
        self.z = z
        self.max1D = max1D

    def move(self):
        self.x += (2 * rr() - 1.0) * self.max1D
        self.y += (2 * rr() - 1.0) * self.max1D
        self.z += (2 * rr() - 1.0) * self.max1D

    def position(self):
        return (self.x, self.y, self.z)


class SphereInteraction:
    __metaclass__ = ABCMeta

    #@abstractmethod
    #def energy(self, group, sphere):
    #    pass

    @abstractmethod
    def energyCallback(self, sphere):
        pass

class LJInteraction(SphereInteraction):
    
    def __init__(self, epsilon, sigma, group, cutoff=10.0):
        self.epsilon = epsilon
        self.sigma = sigma
        self.group = group
        self.cutoff = cutoff
        
    def energyCallback(self, sphere):
        callback = LJPotential(self.epsilon, self.sigma, sphere.position()) 
        return self.group.reduce(callback, sphere.position(), self.cutoff)

    def createCallback(self, sphere):
        # Here a function is created and returned!
        return LJPotential(self.epsilon, self.sigma, sphere.position())


class LJPotential:
    
    def __init__(self, epsilon, sigma, offset):
        self.epsilon = epsilon
        self.sigma = sigma
        self.offset = offset

    def value(self, r):
        s = (self.sigma**2 / ((r[0] - self.offset[0])**2 + \
                                  (r[1] - self.offset[1])**2 + \
                                  (r[2] - self.offset[2])**2))**3
        return 4 * self.epsilon * (s - 1) * s


class Always:
    
    def __init__(self, answer=True):
        self.answer = answer

    def accept(self, e0, e1):
        return self.answer



def main():
    nSolvent = 10
    nSolute = 2
    solventParticles = [Sphere(rr(), rr(), rr()) for i in range(nSolvent)]
    soluteParticles = [Sphere(rr(), rr(), rr()) for i in range(nSolute)]

    for p in solventParticles:
        print p.x, p.y, p.z
    

    group1 = SphereGroup(solventParticles, "solvent")
    group2 = SphereGroup(soluteParticles, "solute")

    epsilon = 1.0
    sigma = 1.0
    lj12 = LJInteraction(epsilon, sigma, group1)
    group2.addInteraction(lj12)
    lj21 = LJInteraction(epsilon, sigma, group2)
    group1.addInteraction(lj21)

    particleGroups = []
    particleGroups.append(group1)
    particleGroups.append(group2)

    group1.move(rule=Always(True))
    group2.move(rule=Always(False))

if __name__ == '__main__':
    main()
