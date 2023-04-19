.. _tablesbasics:

Translate Route Into Tables
===========================

A simple two-step demo synthetic route that will be used in this documentation
is shown below. The utilization of the reagents and solvents is being shown as
molar equivalents and volumes (L/kg), respectively, which would have been
determined in the lab. Prices ($/kg) of some compounds are provided as well.
**In essence, the costcalc2 algorigthm will use this information to
calculate the unknown $/kg (RMC) values for Intermediate A and Product.** As
part of this calculation, process mass intensity (PMI) will also be calculated
as a measure of route waste. Note: All of the information in this figure is
fake and is provided for pedagogical purposes only.   

.. image:: ./_images/br_route.png
   :align: center
   :scale: 25 %

Notes on Pricing
----------------

It is important to understand that the *costcalc* algorithm takes price
information for reagents and solvents and calculates *total raw material costs
(RMCs)* for the intermediates and products. In a simplified sense, these RMCs
are the sum of the raw material prices used to synthesize the route
intermediates/product. However, RMCs are only one component of the final price
of a compound, in addition to labor, facilities charges, waste disposal,
profit margin, etc. So **the RMC for the final product is NOT the expected
sales price** and needs to be interpreted with that knowledge in mind.

Accurate raw material prices are a key component of the calculation, and they
are not always trivial to obtain. Prices are highly dependent on the scale,
region where purchase, date when prices were quoted, etc. As such, no price
information (or material information of any kind) is provided in the
*costcalc* module. 

Tables Overview
---------------

The *costcalc* algorithm, including the `web application
<https://costcalc.rnelsonchem.com/>`_ interface, require the chemical
synthetic route information above to be converted into two different input
tables. A "Materials" table will contain all of the information about the raw
materials, such as molecular weight, density, and price. A "Route" table will
contain the route-specific information, such as equivalents of reagents,
volumes of solvent, etc. These tables require some specific formatting, which
will be outlined here. The web application uses an Excel file as the tabular
input. A :download:`blank template <./_downloads/blank_v2.xlsx>` is provided
here, which contains the proper column names and tabs as well as well as some
conditional formatting to help guide user input. A complete :download:`working
Excel file <./_downloads/working_v3.xlsx>` containing this demo route
information is also available.

Materials Table
---------------

The Materials table will be added to the Materials tab in the Excel file, and
this must be populated first. (This is because of some data validation that is
done to ensure correct inputs in the Route tabs.) For the web application,
there should be only one Materials tab per Excel file and the name of the tab
must remain unchanged. Although this table must contain all of the materials
for the synthetic route to be costed, it can contain any arbitrary number of
additional materials as well. This can be seen in the completed Excel file
that can be downloaded above. 

A snapshot of a completed Materials table for this demo route, along with some
additional materials, is shown below. Detailed column descriptions are given
below, but a few general notes are provide here. First, the order of the
materials in this table is not important. In this demo, the compounds with
known prices were separated from compounds without known prices for clarity.
Second, blank lines and formatting are ignored, which can be helpful to make
the input file a little easier to understand. In this case, a blank line
between known/unknown compounds has been added, and clarifying comments have
been added to the column name cells. However; be careful about adding
inadvertent information into "blank" lines. For example, if information was
added into cell E9 below, the costing calculation will fail with an error.
Some of these notes do not apply to the Python programming interface to the
*costcalc* algorithm, as will be discussed in later sections.   

.. image:: ./_images/materials.png
   :align: center


The Materials table contains the following columns, which must be present.
Additional user-defined columns can be added, but will be ignored by the
*costcalc* algorithm.

* *Compound*: This is the desired name for the compounds in the route. These
  names must be unique, i.e. there can not be two compounds with the same
  name.

* *MW*: The molecular weight of the compound in kg/kmol (equivalent to
  g/mol). Every compound must have a MW, so if a compound name is given, the
  MW column will be highlighted automatically to indicate that this value must
  be provided. This value can be a little tricky for materials like celite,
  which may not have a defined MW. Some notes for corner cases like this is
  provided below.

* *Density*: The density of the compound in kg/L (equivalent to g/mL). This
  value is only mandatory for compounds that will be used as solvents, as
  described in the Route table section, otherwise they can be left blank. If a
  density is added for a non-solvent reagent, it will be ignored.

* *$/kg*: The price of the compound in $/kg. Most compounds will need to have
  a defined price with the following exception. If the compound represents a
  product of a reaction, then the price is optional. In these cases a price
  can be provided, if desired, but it may be ignored during the calculations. 

* *Notes*: (Optional) Notes for the compound information. For example, this is
  a good place to add the date and source of the pricing information and the
  name of the person who acquired that information.

Route Table
-----------

The route information is added into a separate tab in the Excel file with the
required columns described in the list below. There can be one or more route
tabs per Excel file, and the tab names can be arbitrarily chosen (except for
the name "Materials" as described above). In the blank template file, the
route tab is called "Route 1"; in the working demo file, the route tab has
been renamed to "Bromine Route". If a new route tab is desired, the column
names are the most critical component that needs to be copied. However, the
route tabs in the Excel files provided here do contain some conditional
formatting and data validation to ensure that values are added correctly, so
it is recommended to use the "Move or copy..." dialog (accessed by
right-clicking on the tab name) to duplicate a current route tab.

A snapshot of the Route table for the demo route above is shown below.
Pictures (ChemDraw or regular images files) of the route can be added to this
tab, and they will be ignored by the *costcalc* algorithm. Descriptions of the
columns are provided in the list below the figure. As with the Materials
table, additional user-defined columns or blank lines can be added, but they
will be ignored during the costing operation.

.. image:: ./_images/reactions.png
   :align: center

* *Step*: An unique identifier to delineate the synthetic step in the route.
  These can be simply numerical numbers (e.g. 1, 2, 3) and/or text ("1a" or
  "Int A"). Steps do not need to be added into the table in any particular
  order, as they will be automatically sorted during the costing calculation.
  In fact, the compounds from every step could be added in arbitrary order;
  however, this is not recommended from a clarity standpoint.

* *Compound*: The name of a reagent/solvent/product for the step. These names
  must *exactly* correspond to the Materials table, so a drop-down selector is
  provided to ensure that a valid name is selected. (This is why the Materials
  table should be created first.)

* *Equiv*: Molar equivalents of a reagent or product. Although this value can
  be used for solvents, it is more common to define solvent utilization with
  *Volumes*, as described in the next column. These values can be scaled as
  needed, but they are typically scaled such that the limiting reagent is 1
  equivalent. For a product, the equivalents are the theoretical equivalents
  multiplied by the fractional percent yield. For example, in a reaction with
  a starting material to product ratio of 1:1 and a 75% yield of product, the
  equivalents of product would be :math:`1*0.75=0.75`. If 2 moles of product
  are expected (e.g. breaking up a dimer) with the same reaction yield, the
  equivalents of product would be :math:`2*0.75=1.5`.

* *Volumes*: The amount of solvent utilization in volumes. This value is only
  required if *Equiv* for a particular compound is not given; if this column
  is used, the next two columns (*Relative* and *Sol Recyc*) are required. The
  unit for volumes is L/kg, which can be interpreted as "liters of this
  solvent per kg of a reference compound." This is numerically equivalent to
  mL/g. The reference compound is defined in the next column. 

* *Relative*: The reference compound for solvent volume calculations. This is
  typically the starting material/limiting reagent of the reaction, but that
  may not always be the case. Again, the name here must correspond to a
  compound from the Materials table; this cell contains a drop-down selector to
  ensure that a valid compound name is selected. The material name must also
  be defined in the current reaction *Step*, otherwise the cost calculation
  will result in an error.

* *Sol Recyc*: The fractional percentage of this solvent that it is expected
  could be recycled. For example, if 95% of the solvent can be recycled, then
  this cell will contain the value 0.95. In our demo example, we are assuming
  that 75% of the solvents can be recycled; however, if you are unsure, set
  this value to 0, which means that none (0%) of this solvent can be recycled.

* *Cost Step*: The step identifier that indicates where the RMC for this
  compound will be calculated. The value here must be a valid entry from the
  *Step* column, and these entries are only necessary for route intermediates
  and the overall product. (I.e. any compound that does not have a $/kg entry
  in the Materials table.) This column is critical as it provides a "roadmap"
  of sorts to define how the different reactions are connected. In our demo
  example, the RMC for Intermediate A is calculated in step "1", so *all
  usages of Intermediate A must be labeled as "1"*. A simplified version of a
  longer linear and convergent route are provide below for additional
  demonstration purposes.

.. _OPEXinput:

* *OPEX*: (Optional) An estimate, in $/kg, of the operating expenses for a
  given reaction step. This number is only valid for the product of any given
  step.  Although these values are not given for the current demo route, they
  could have been given for Intermediate A in Step 1 (Cell H5) and/or Product in
  Step 2 (cell H10). For route intermediates, these values are added to the
  RMC values in subsequent steps; the OPEX for the final product is added to
  the final RMC value in the *$/kg* column. This can be a bit confusing at
  first, so a :ref:`second model using OPEX values <OPEX>` will be presented
  in the next section. 

* *Notes*: (Optional) Notes for this particular compound. For example, a
  reference can be included here if the reaction was taken from the
  literature, or a short bit of text can be added to acknowledge any
  assumptions in the numbers.

Linear vs Convergent Syntheses
______________________________

The *Step* and *Cost Step* columns and their connections are vital to ensure
that the route is costed correctly. Using these connection schemes we can also
define routes of arbitrary number of steps and level of convergence. Below are
two very simplified Route tables for the products of two different three-step
synthetic routes. One is completely linear and the other is convergent. 

Below is the simplified Route table for the products of a three step linear
route, which is shown in the figure as well. The identifiers in the *Step* and
*Cost Step* columns have been color coded for additional clarity.

.. image:: ./_images/3step_linear.png
   :align: center

The next figure is a simplified Route table for the products of a convergent
three-step route, as shown. Again, color coding is added to clarify values
that must be the same.

.. image:: ./_images/3step_conv.png
   :align: center




