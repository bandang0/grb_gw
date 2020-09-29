'''No idea.'''

from helpers import *

print(polyfit(log(OBS_TIMES_1 / 100.), log(OBS_FLUXES_1), OBS_ERRORS_1 / OBS_FLUXES_1))
print(polyfit(log(OBS_TIMES_2 / 300.), log(OBS_FLUXES_2), OBS_ERRORS_2 / OBS_FLUXES_2))
