Tutorial
++++++++

Creating a new sample
---------------------
A sample is the basic structure of the package. In order to analyze your data you will have to create a new sample first.
lets start with importing the package.

.. code-block:: python

   from RockPyV3 import Sample

This will import the needed parts of the module, so they can be used.

Now lets create a new sample.

.. code-block:: python

   first_sample = Sample(name='sample 01')

The sample class has a lot more to offer. You can look it up here :doc:`structure`. But let's add some more information,
 like the mass to the sample, so we can actually normalize the data properly.

.. code-block:: python

   first_sample = Sample(name='sample 01', mass=35.2, mass_unit='mg')

As you can see we added the mass and its unit. The mass is stored in kg. If no unit is added, the software assumes you
entered **mg**.

If you want to change or add some more information, you can do it by using the following methods:

- first_sample.add_mass(mass, mass_unit)
- first_sample.add_height(height, length_unit)
- first_sample.add_diameter(diameter, length_unit)

.. warning:: length_unit is global, if you specify a different length_unit for **height** and **diameter** you may get into trouble


adding a measurement
====================

So we created the first sample. Since we have already measured a hysteresis loop, we can add this measurement to the sample.

.. code-block:: python

   first_sample.add_measurement(mtype, mfile, machine)

:code:`mtype` stands for **measurement type**. In this case it would be a hysteresis, therefore :code:`mtype='hys'`

:code:`mfile` stands for **measurement file** and is the location of the file on your harddrive.

:code:`machine` is the **machine** the measurement was done on. We will use the VSM/AGM from PMC in this example.

::

   hysteresis = first_sample.add_measurement(mtype='hys', mfile='\User\data\first_sample.hys', machine='vsm')

So now we not only added a hysteresis to the sample, we also created a hysteresis object, which we can now use to determine the usual hysteresis parameters. More on this later or here :doc:`measurements`

adding a treatment to a measurement
===================================

Creating a sample group
-----------------------
