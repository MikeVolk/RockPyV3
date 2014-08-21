Measurement
-----------

The measurement class is the parent class for all child measurements (e.g. Hysteresis). Its constructor needs
a sample object (sample_obj), a measurement type (mtype), measurement file (mfile) and the machine (machine) the measurement has been
done with. Optionally a name for the logger can be passed and will be initialized, so the logging shows
'RockPy.MEasurement.Hysteresis' instead of 'RockPy.MEASURENT'.

Derived classes need to initialize the Measurement class with:

:code:

    super(Derived_Class, self).__init__(sample_obj, mtype, mfile, machine, log)

but may have aditional constructor parameters. Example:

:code:

   class Derived_Class(Measurement):
           def __init__(self, sample_obj,
                        af_obj, parm_obj,
                        mtype=None, mfile=None, machine=None,
                        **options):

In this example the parent constructor does not get the af_obj or the parm_obj. The **option is needed
for dynamic construction of the class.

.. autoclass:: Structure.measurements.Measurement
:members:

Hysteresis
++++++++++

.. autoclass:: Structure.measurements.Hysteresis
:members:

Paleointensities
++++++++++++++++

.. autoclass:: Structure.measurements.Thellier
:members: