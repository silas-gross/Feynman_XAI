from sympy import *
free_1, free_2, free_3, free_4=symbols('free_1 free_2 free_3 free_4')

#integrate(1/(free_1**2 + 156.25)*1/(free_2**2 + 15600), (free_1, 0, 637.2644), (free_2, 0, 637.2644), (free_3, 0, 637.2644), (free_4, 0, 637.2644))
#integrate(1/(free_1**2 + 156)*1/(free_2**2 + 15775)*1/(free_3**2 + 156), (free_1, 0, 637.2644), (free_2, 0, 637.2644), (free_3, 0, 637.2644), (free_4, 0, 637.2644))

integrate(1/(free_1**2 + 156)*1/(free_2**2 + 15775)*1/(free_3**2 + 156)*1/(free_4**2 + 15775), (free_1, 0, 637.2644), (free_2, 0, 637.2644), (free_3, 0, 637.2644), (free_4, 0, 637.2644))

