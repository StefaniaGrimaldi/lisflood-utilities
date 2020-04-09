import os
import logging
import itertools

from nine import IS_PYTHON2
if IS_PYTHON2:
    from pathlib2 import Path
else:
    from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from lisfloodutilities.readers import PCRasterMap


logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger()


class Comparator(object):
    """

    """
    glob_expr = None

    def __init__(self, array_equal=False, for_testing=True):
        """

        """
        self.array_equal = array_equal
        self.for_testing = for_testing
        self.errors = []

    def compare_dirs(self, path_a, path_b, skip_missing=True):
        logger.info('Comparing %s and %s [skip missing: %s]', path_a, path_b, str(skip_missing))
        path_a = Path(path_a)
        path_b = Path(path_b)
        for fa in itertools.chain(*(path_a.glob(e) for e in self.glob_expr)):
            fb = path_b.joinpath(fa.name)
            if not fb.exists():
                if skip_missing:
                    logger.info('skipping %s as it is not in %s', fb.name, path_b.as_posix())
                    continue
                else:
                    message = '{} is missing in {}'.format(fb.name, path_b.as_posix())
                    if self.for_testing:
                        assert False, message
                    else:
                        self.errors.append(message)

            self.compare_files(fa.as_posix(), fb.as_posix())
        if not self.for_testing:
            return self.errors

    def compare_files(self, fa, fb):
        raise NotImplementedError()


class PCRComparator(Comparator):
    glob_expr = ['**/*.[0-9][0-9][0-9]', '**/*.map']

    def compare_files(self, file_a, file_b):
        logger.info('Comparing %s and %s', file_a, file_b)
        map_a = PCRasterMap(file_a)
        map_b = PCRasterMap(file_b)
        if map_a != map_b:
            map_a.close()
            map_b.close()
            message = '{} different from {}'.format(file_a, file_b)
            if self.for_testing:
                assert False, message
            else:
                self.errors.append(message)
                return message
        else:
            map_a.close()
            map_b.close()
            assert True


class TSSComparator(Comparator):
    glob_expr = ['**/*.tss']

    def compare_files(self, file_a, file_b):
        logger.info('Comparing %s and %s', file_a, file_b)
        with open(file_a, 'rb') as fp1, open(file_b, 'rb') as fp2:
            # need to skip first line in TSS as it reports settings filename
            fp1.readline()
            fp2.readline()
            while True:
                b1 = fp1.readline()
                b2 = fp2.readline()
                if (b1 != b2) or (not b1 and b2) or (not b2 and b1):
                    message = '{} different from {}'.format(file_a, file_b)
                    if self.for_testing:
                        assert False, message
                    else:
                        self.errors.append(message)
                        return message
                if not b1:
                    return None


class NetCDFComparator(Comparator):
    glob_expr = ['**/*.nc']

    def __init__(self, mask, atol=0.05, rtol=0.1, max_perc_diff=0.2, max_perc_large_diff=0.1,
                 array_equal=False, for_testing=True):
        super(NetCDFComparator, self).__init__(array_equal=array_equal, for_testing=for_testing)
        if isinstance(mask, str):
            mask = Dataset(mask)
            maskvar = [k for k in mask.variables if len(mask.variables[k].dimensions) == 2][0]
            self.maskarea = np.logical_not(mask.variables[maskvar][:, :])
        else:
            self.maskarea = mask
        self.atol = atol
        self.rtol = rtol
        self.max_perc_large_diff = max_perc_large_diff
        self.max_perc_diff = max_perc_diff
        self.large_diff_th = self.atol * 10

    def compare_arrays(self, vara_step, varb_step, varname=None, step=None, filepath=None):

        vara_step = np.ma.compressed(np.ma.masked_array(vara_step, self.maskarea)).astype('float64')
        varb_step = np.ma.compressed(np.ma.masked_array(varb_step, self.maskarea)).astype('float64')
        diff_values = np.ma.abs(vara_step - varb_step)
        diff_values = diff_values[~np.isnan(diff_values)]
        same_values = np.ma.allclose(diff_values, np.zeros(diff_values.shape), atol=self.atol, rtol=self.rtol)
        all_ok = vara_step.size == varb_step.size and same_values
        array_ok = np.isclose(diff_values, np.zeros(diff_values.shape), atol=self.atol, rtol=self.rtol, equal_nan=True)
        different_values_size = array_ok[~array_ok].size
        if (not all_ok) and (different_values_size > 0):
            max_diff = np.ma.amax(diff_values)  # returns a scalar
            perc_wrong = different_values_size * 100 / vara_step.size
            result = np.ma.where(diff_values >= max_diff)
            rel_diff = max_diff * 100. / np.maximum(vara_step[result], varb_step[result])
            if rel_diff.size > 0 and np.max(rel_diff) > 0.01 and (perc_wrong >= self.max_perc_diff or (perc_wrong >= self.max_perc_large_diff and max_diff > self.large_diff_th)):
                step = step if step is not None else '(no time)'
                filepath = os.path.basename(filepath) if filepath else '<mem>'
                varname = varname or '<unknown var>'
                message = '{}/{}@{} - {:3.2f}% of different values - max diff: {:3.2f}'.format(filepath, varname, step, perc_wrong, max_diff)
                logger.error(message)
                if self.for_testing:
                    assert False, message
                else:
                    self.errors.append(message)
                    return message
        assert True

    def compare_files(self, file_a, file_b):
        logger.info('Comparing %s and %s', file_a, file_b)
        with Dataset(file_a) as nca, Dataset(file_b) as ncb:
            num_dims = 3 if 'time' in nca.variables else 2
            var_name = [k for k in nca.variables if len(nca.variables[k].dimensions) == num_dims][0]
            vara = nca.variables[var_name]
            try:
                varb = ncb.variables[var_name]
            except KeyError:
                var_nameb = [k for k in ncb.variables if len(ncb.variables[k].dimensions) == num_dims][0]
                message = 'Files: {} vs {} have different variables names A:{} B:{}'.format(file_a, file_b, var_name, var_nameb)
                if self.for_testing:
                    assert False, message
                else:
                    self.errors.append(message)
                    return message
            if 'time' in nca.variables:
                stepsa = nca.variables['time'][:]
                stepsb = ncb.variables['time'][:]
                if len(stepsa) != len(stepsb):
                    len_stepsa = len(stepsa)
                    len_stepsb = len(stepsb)
                    message = 'Files: {} vs {}: different size in time axis A:{} B:{}'.format(file_a, file_b, len_stepsa, len_stepsb)
                    if self.for_testing:
                        assert False, message
                    else:
                        self.errors.append(message)
                        return message
                for step, _ in enumerate(stepsa):
                    values_a = vara[step][:, :]
                    values_b = varb[step][:, :]
                    if not self.array_equal:
                        mess = self.compare_arrays(values_a, values_b, var_name, step, file_a)
                    else:
                        mess = self.compare_arrays_equal(values_a, values_b, var_name, step, file_a)
                    if mess:
                        self.errors.append(mess)
            else:
                values_a = vara[:, :]
                values_b = varb[:, :]
                if not self.array_equal:
                    mess = self.compare_arrays(values_a, values_b, var_name, filepath=file_a)
                else:
                    mess = self.compare_arrays_equal(values_a, values_b, var_name, filepath=file_a)
                if mess:
                    self.errors.append(mess)
            return self.errors

    def compare_arrays_equal(self, values_a, values_b, varname=None, step=None, filepath=None):
        filepath = os.path.basename(filepath) if filepath else '<mem>'
        message = '{}/{}@{} is different'.format(filepath, varname, step)
        if self.for_testing:
            np.testing.assert_array_equal(values_a, values_b, message)
        else:
            if not np.array_equal(values_a, values_b):
                self.errors.append(message)
                return message
