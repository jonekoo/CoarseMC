module m_rod_interaction_inc
  use m_rod, t => rod, u => rod
  use m_particlegroup, group => particlegroup
  use m_stopping_interaction
  implicit none

  include 'tu_interaction-def.inc'

contains

  include 'tu_interaction-proc.inc'
  
end module m_rod_interaction_inc


module m_rod_interaction
  use m_rod_interaction_inc, rod => t, rod => u, &
       rod_interaction => tu_interaction, rodrod_potential => tu_potential
end module m_rod_interaction
