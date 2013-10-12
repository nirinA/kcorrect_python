import kcorrect
import numpy
import os
test_file = os.path.join(os.environ['KCORRECT_DIR'], 'test', 'sample.dat')
kcorrect.load_templates()
kcorrect.load_filters()
kcorrect.fit_coeffs_from_file(test_file)
kcorrect.reconstruct_maggies_from_file('coeffs.dat')
