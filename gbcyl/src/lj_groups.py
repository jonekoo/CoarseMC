#-*- coding: utf-8 -*-
from abc import ABCMeta, abstractmethod

class ParticleGroup:
    __metaclass__ = ABCMeta

    def __init__(self, particles, id):
        self.particles = particles

    @abstractmethod
    def move(self):
        pass

    @abstractmethod:
    def add(particle):
        pass


class SphereGroup(ParticleGroup):
    
    def move(self, rule):
        for p in self.particles:
            p0 = copy(p)
            e0 = self.energy(p)
            p.move()
            e1 = self.energy(p)
            if rule.accept(e0, e1):
        return 

    def energy(self, particle):
        for ia in self.interactions:
            energy = energy + ia.energy(particle)
        return energy

    def getNeighbours(position, cutoff):
        nbs = []
        for p in self.particles:
            if abs(p.position - position) < cutoff:
                nbs.append(p)
        return nbs

    def potentialThroughNeighbours(position, cutoff, potential):
        sum = 0
        for p in self.particles:
            if abs(p.position - position) < cutoff:
                # Jälkimmäinen argumentti callback-kutsussa ehkä turha
                sum = sum + potential.value(p.position)
        return sum


class SphereInteraction:
    __metaclass__ = ABCMeta

    @abstractmethod
    def energy(self, group, sphere):
        pass


class LJInteraction(SphereInteraction):
    
    def __init__(self, epsilon, sigma, group, cutoff=None):
        self.epsilon = epsilon
        self.sigma = sigma
        self.group = group
        
    def energy(self, sphere):
        neighbours = self.group.getNeighbours(sphere.position, self.cutoff)
        for nb in neighbours:
            energy = energy + self.potential(epsilon, sigma, nb.position - sphere.position)

    def energyCallback(self, sphere):
        callback = self.createCallback(sphere)
        group.callbackThroughNeighbours(sphere.position, self.cutoff, callback)

    def potential(epsilon, sigma, r):
        s = (sigma / r)**6
        return 4 * epsilon * (s - 1) * s

    def createCallback(self, sphere):
        # Here a function is created and returned!
        return LJPotential(self.epsilon, self.sigma, sphere.position)


def LJPotential:
    
    def __init__(epsilon, sigma, offset):
        self.epsilon = epsilon
        self.sigma = sigma
        self.offset = offset

    def value(self, r):
        s = (sigma / (r - offset))**6
        return 4 * epsilon * (s - 1) * s

