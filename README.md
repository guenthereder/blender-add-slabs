# blender-add-slabs

This blender add-on creates the menu item "add slabs" under object. Selected a polygon "add slabs" creates half-plane slabs with an inclination of 45 degree that originate at the polygon edges. These slabs are bounded by the bisectors of the polygon edges (also lifted to the same inclination). This add-on helps visualize certain straight skeleton approaches where such a specific half-plane arrangement is used.

## Weights

Multiplicative weights are considered if they are part of the input .obj file. A weight for an edge (v1,v2) can be given as an attribute for v1. 

Example:

v x-value y-value z-value
v x-value y-value z-value
v x-value y-value z-value
vt -100 0.5
vt -100 0.7
f 1/2 2/1 3/1

In this example we give three vertices (v) and define a triangle (with f). In the listing of the vertices used for this face (f) we also reference the used attribute. Then the first edge of vertex 1 to 2 has attribute 2, the other edges have attribute 1. The first attribute value is not relevant, the weight is given with the last value. The weight range goes from 0 to 1.
