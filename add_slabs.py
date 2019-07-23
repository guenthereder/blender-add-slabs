"""
    add_slabs.py is a Blender Add-On that computes and constructs edge and
    'motorcycle' slabs from a given input polygon.
    {one line to give the program's name and a brief idea of what it does.}
    Copyright (C) 2016  Günther Eder  geder@cosy.sbg.ac.at

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


bl_info = {
    "name": "Add Edge/MC Slabs",
    "author": "Günther Eder",
    "version": (2, 1),
    "blender": (2, 80, 0),
    "api": 33333,
    "location": "Object > Add Slabs",
    "description": "Adding Edge (and Motorcylce) Slabs",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Add Mesh"}

import bpy, math, random, bmesh
import numpy as np
from mathutils import Vector, Euler, Matrix

class AddEdgeSlabs(bpy.types.Operator):
    """Add Edge/MC Slabs"""
    bl_idname = "object.addslabs"
    bl_label = "Add Slabs"
    bl_options = {'REGISTER', 'UNDO'}

    total = bpy.props.IntProperty(name="Slope (Deg)", default=45, min=1, max=89)

    MAX_ROTATION = math.pi/2.0
    maxEdgeLength = 0.0

    def unit_vector(self, vector):
        """ Returns the unit vector of the vector.  """
        return vector / np.linalg.norm(vector)

    def perp( self, a ) :
        b = np.empty_like(a)
        b[0] = -a[1]
        b[1] = a[0]
        return b

    """ http://stackoverflow.com/questions/3252194/numpy-and-line-intersections """
    # line segment a given by endpoints a1, a2
    # line segment b given by endpoints b1, b2
    def seg_intersect(self, a1,a2, b1,b2) :
        da = a2-a1
        db = b2-b1
        dp = a1-b1
        dap = self.perp(da)
        denom = np.dot( dap, db)
        num = np.dot( dap, dp )
        return (num / denom.astype(float))*db + b1



    """ http://stackoverflow.com/questions/31735499/calculate-angle-clockwise-between-two-points """
    def angle_between(self, v1, v2):
        p1 = [v1.x,v1.y]
        p2 = [v2.x,v2.y]
        ang1 = np.arctan2(*p1[::-1])
        ang2 = np.arctan2(*p2[::-1])
        return (ang1 - ang2) % (2 * np.pi)

    def get_random_color(self):
        ''' generate rgb using a list comprehension '''
        r, g, b = [random.random() for i in range(3)]
        return r, g, b, 1

    """ 
    Rand Material :
        https://blenderartists.org/forum/showthread.php?210937-How-do-I-add-a-new-material-using-python
    """ 
    def returnMaterialByName(self, passedName):
        result = None
        for m in bpy.data.materials:
            if m.name == passedName:
                result = m
                break
        return result

    def createNewMaterial (self, passedName):
        tempMat = bpy.data.materials.new(passedName)
        if tempMat != None:
            r = random.random()
            g = random.random()
            b = random.random()
            a = 1.0
            tempMat.diffuse_color = (r,g,b,a)
            #tempMat.diffuse_shader = 'LAMBERT'
            #tempMat.diffuse_intensity = 1.0
            tempMat.specular_color = (0.9,0.9,0.9)
            #tempMat.specular_shader = 'COOKTORR'
            #tempMat.specular_intensity = 0.5
            #tempMat.alpha = 1.0
            #tempMat.ambient = 0.3
            #tempMat.emit = 0.2
        return tempMat

    def aquireOrCreateMaterial(self, passedName):
        tempMat = self.returnMaterialByName(passedName)
        if tempMat == None:
            tempMat = self.createNewMaterial(passedName)
        return tempMat

    def fetchIfObject (self, passedName= ""):
        try:
            result = bpy.data.objects[passedName]
        except:
            result = None
        return result

    def rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        source: unutbu in http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
        """
        axis = np.asarray(axis)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


    def add_slab_geometry(self, context,planeVertices,planeFaces,cnt):
        meshName     = "Plane %s" % (cnt)
        if len(planeVertices) < 4:
            meshName += " MC"
        materialName = "Material (rand) %s" % (cnt)
        mymesh   = bpy.data.meshes.new(meshName)
        myobject = bpy.data.objects.new(meshName, mymesh)
        myobject.color = self.get_random_color()

        #myobject.location = bpy.context.scene.cursor_location
        #context.scene.objects.link(myobject)
        context.collection.objects.link(myobject)

        #Create mesh
        mymesh.from_pydata(planeVertices,[],planeFaces)

        me  = myobject.data
        mat = self.aquireOrCreateMaterial(materialName)
        me.materials.append(mat)
                
        mymesh.update(calc_edges=True)
        
        return

    def add_slab(self, context,v0,v1,v2,v3,cnt,rotation):
        # direction vector from 0 -> 1 and 1 -> 0
        aLeftDir  = self.unit_vector(v0 - v1)
        aRightDir = self.unit_vector(v2 - v1)
        bLeftDir  = self.unit_vector(v1 - v2)
        bRightDir = self.unit_vector(v3 - v2)

        """ angle is calculated clockwise """
        aDir = self.angle_between(aLeftDir,aRightDir)
        bDir = self.angle_between(bLeftDir,bRightDir)

        eDir = self.unit_vector(v2 - v1)

        # rotate them by 90 deg (2D) + norm and unit
        eNewDir = Vector((-eDir.y, eDir.x, 0))

        rotation_axis = (v1 - v2)
        rot_mat = self.rotation_matrix(rotation_axis,rotation)
        eNewDir = Vector(np.dot(rot_mat, eNewDir))
        eNewDir *= self.maxEdgeLength
        
        # two extra vertices for plane
        v4 = (v2 + eNewDir)
        v5 = (v1 + eNewDir)

        planeVertices = [(v1.x,v1.y,v1.z),(v2.x,v2.y,v2.z),(v4.x,v4.y,v4.z),(v5.x,v5.y,v5.z)]
        planeFaces = [(0,1,2,3)]
        self.add_slab_geometry(context,planeVertices,planeFaces,cnt)

        if aDir > math.pi:
            aRotAxis = Vector((0,0,1))
            aRotMat  = self.rotation_matrix(aRotAxis,aDir/2.0)
            mcDir    = Vector(np.dot(aRotMat, aRightDir))
            
            alpha = (math.pi-aDir)/2.0
            mcLen = self.maxEdgeLength/math.cos(alpha)

            baseLen = math.fabs(mcLen * math.sin(alpha))
        
            vm = (v5 - (aRightDir * baseLen * 2))
            v_on_mc = (v1 + (mcDir * self.maxEdgeLength * 2))

            p1 = np.array( [vm.x,vm.y] )
            p2 = np.array( [v5.x,v5.y] )

            p3 = np.array( [v1.x,v1.y] )
            p4 = np.array( [v_on_mc.x,v_on_mc.y] )

            vm_xy = self.seg_intersect( p1,p2, p3,p4)
        
            vm = Vector((vm_xy[0],vm_xy[1],vm.z))

            planeVerticesM = [(v1.x,v1.y,v1.z),(v5.x,v5.y,v5.z),(vm.x,vm.y,vm.z)]
            planeFacesM = [(0,1,2)]
            self.add_slab_geometry(context,planeVerticesM,planeFacesM,cnt)


        if bDir > math.pi:
            bRotAxis = Vector((0,0,-1))
            bRotMat  = self.rotation_matrix(bRotAxis,bDir/2.0)
            mcDir    = Vector(np.dot(bRotMat, bLeftDir))

            alpha = (math.pi-bDir)/2.0
            mcLen = self.maxEdgeLength/math.cos(alpha)

            baseLen = math.fabs(mcLen * math.sin(alpha))

            vm = (v4 - (bLeftDir * baseLen * 2))
            v_on_mc = (v2 + mcDir * self.maxEdgeLength * 2)

            p1 = np.array( [vm.x,vm.y] )
            p2 = np.array( [v4.x,v4.y] )

            p3 = np.array( [v2.x,v2.y] )
            p4 = np.array( [v_on_mc.x,v_on_mc.y] )

            vm_xy = self.seg_intersect( p1,p2, p3,p4)
        
            vm = Vector((vm_xy[0],vm_xy[1],vm.z))
            
            planeVerticesM = [(v2.x,v2.y,v2.z),(vm.x,vm.y,vm.z),(v4.x,v4.y,v4.z)]
            planeFacesM = [(0,1,2)]
            self.add_slab_geometry(context,planeVerticesM,planeFacesM,cnt)
        
        return


    def getMaxEdgeLengthAndBB(self, Obj, OWMatrix):
        maxEdgeLength = 0
        min_x = float("inf")
        min_y = float("inf")
        max_x = -float("inf")
        max_y = -float("inf")
        # find the maximal edge Length in input
        for e in Obj.data.edges :
            v0_idx = e.vertices[0]
            v1_idx = e.vertices[1]
            v0 = OWMatrix @ Obj.data.vertices[v0_idx].co
            v1 = OWMatrix @ Obj.data.vertices[v1_idx].co
            edgeLength = (v0 - v1).length

            if edgeLength > maxEdgeLength:
                maxEdgeLength = edgeLength

            if v0.x < min_x:
                min_x = v0.x
            
            if v0.x > max_x:
                max_x = v0.x
                
            if v0.y < min_y:
                min_y = v0.y
               
            if v0.y > max_y:
                max_y = v0.y

        return maxEdgeLength, min_x, min_y, max_x, max_y

    def addComplementaryPolygon(self, vertices,faces):
        meshName = "Complementary Polygon"
        materialName = "CP Material (rand)"
        mymesh   = bpy.data.meshes.new(meshName)
        myobject = bpy.data.objects.new(meshName, mymesh)
        myobject.color = self.get_random_color()
                
        #Set location and scene of object
        #bpy.context.scene.objects.link(myobject)
        bpy.context.collection.objects.link(myobject)
     
        #Create mesh
        mymesh.from_pydata(vertices,[],faces)

        me  = myobject.data
        mat = self.aquireOrCreateMaterial(materialName)
        me.materials.append(mat)

        mymesh.update(calc_edges=True)

        return


    def execute(self, context):
        SO = context.selected_objects
        Obj = context.active_object
        OWMatrix = Obj.matrix_world

        mesh = Obj.data
        vertices = []
        edgeCnt = 0

        if SO :
            self.maxEdgeLength, min_x, min_y, max_x, max_y = self.getMaxEdgeLengthAndBB(Obj,OWMatrix)


            # CW arranged!
            BBox = [(min_x-(self.maxEdgeLength),min_y-(self.maxEdgeLength),0),
                    (min_x-(self.maxEdgeLength),max_y+(self.maxEdgeLength),0),
                    (max_x+(self.maxEdgeLength),max_y+(self.maxEdgeLength),0),
                    (max_x+(self.maxEdgeLength),min_y-(self.maxEdgeLength),0)]

            ################# add edge slabs ################################
            for face in mesh.polygons:
                for vert in face.vertices:
                    vertices.append(vert)
                    
                    if len(vertices) < 4:
                        continue
                
                    v0 = OWMatrix @ Obj.data.vertices[vertices[-4]].co
                    v1 = OWMatrix @ Obj.data.vertices[vertices[-3]].co
                    v2 = OWMatrix @ Obj.data.vertices[vertices[-2]].co
                    v3 = OWMatrix @ Obj.data.vertices[vertices[-1]].co
                    
                    """ we use the uv y-value of the first vertex that defines 
                        the current slab, e.g., for slab v1v2 we us the uv y-val
                        of v1 
                    """
                    v1_idx = vertices[-3]
                    uv_weight = Vector((-100,0.5))
                    for uv_layer in Obj.data.uv_layers:
                        uv_weight = uv_layer.data[v1_idx].uv    
                    weight = math.pi/2.0 - (uv_weight.y * self.MAX_ROTATION)
                
                    edgeCnt = edgeCnt + 1
                    self.add_slab(context,v0,v1,v2,v3,edgeCnt,weight)
                    
                    
                #first edges:
                v0 = OWMatrix @ Obj.data.vertices[vertices[-1]].co
                v1 = OWMatrix @ Obj.data.vertices[vertices[0] ].co
                v2 = OWMatrix @ Obj.data.vertices[vertices[1] ].co
                v3 = OWMatrix @ Obj.data.vertices[vertices[2] ].co
            
                v1_idx = vertices[0]
                uv_weight = Vector((-100,0.5))
                for uv_layer in Obj.data.uv_layers:
                    uv_weight = uv_layer.data[v1_idx].uv    
                weight = math.pi/2.0 - (uv_weight.y * self.MAX_ROTATION)
            
                edgeCnt = edgeCnt + 1
                self.add_slab(context,v0,v1,v2,v3,edgeCnt,weight)
            
                #last edge-1
                v0 = OWMatrix @ Obj.data.vertices[vertices[-3]].co
                v1 = OWMatrix @ Obj.data.vertices[vertices[-2]].co
                v2 = OWMatrix @ Obj.data.vertices[vertices[-1]].co
                v3 = OWMatrix @ Obj.data.vertices[vertices[0] ].co
            
                v1_idx = vertices[-2]
                uv_weight = Vector((-100,0.5))
                for uv_layer in Obj.data.uv_layers:
                    uv_weight = uv_layer.data[v1_idx].uv    
                weight = math.pi/2.0 - (uv_weight.y * self.MAX_ROTATION)
            
                edgeCnt = edgeCnt + 1
                self.add_slab(context,v0,v1,v2,v3,edgeCnt,weight)
                
                #last edge
                v0 = OWMatrix @ Obj.data.vertices[vertices[-2]].co
                v1 = OWMatrix @ Obj.data.vertices[vertices[-1]].co
                v2 = OWMatrix @ Obj.data.vertices[vertices[0] ].co
                v3 = OWMatrix @ Obj.data.vertices[vertices[1] ].co
            
                v1_idx = vertices[-1]
                uv_weight = Vector((-100,0.5))
                for uv_layer in Obj.data.uv_layers:
                    uv_weight = uv_layer.data[v1_idx].uv    
                weight = math.pi/2.0 - (uv_weight.y * self.MAX_ROTATION)
            
                edgeCnt = edgeCnt + 1
                self.add_slab(context,v0,v1,v2,v3,edgeCnt,weight)
            
            ###################################################################

            # add complement of polygon to visually remove slabs outside of P
            cPolyVertices = []
            
            # CW Arranged
            idx = 4
            faceList = []
            bboxAdded = False

            # BBox vertices
            for v in BBox:
                cPolyVertices.append(v)

            for face in mesh.polygons:
                for v_idx in face.vertices:
                    v = OWMatrix @ Obj.data.vertices[v_idx].co
                    cPolyVertices.append(v)
                    faceList.append(idx)
                    if v.x == min_x and (v.y == min_y or v.y == max_y) and not bboxAdded: 
                        if v.y == min_y:
                            faceList.extend([0,1,2,3,0])
                        else:
                            faceList.extend([1,2,3,0,1])

                        faceList.append(idx)
                        bboxAdded = True
                    
                    idx = idx + 1

            faceList = [(faceList)]
            self.addComplementaryPolygon(cPolyVertices,faceList)

        else:
            print("Select a polygon to use this script!")
       
        return {'FINISHED'}


def menu_func(self, context):
    self.layout.operator(AddEdgeSlabs.bl_idname)

# store keymaps here to access after registration
addon_keymaps = []


def register():
    bpy.utils.register_class(AddEdgeSlabs)
    bpy.types.VIEW3D_MT_object.append(menu_func)

    # handle the keymap
    wm = bpy.context.window_manager
    km = wm.keyconfigs.addon.keymaps.new(name='Object Mode', space_type='EMPTY')
    kmi = km.keymap_items.new(AddEdgeSlabs.bl_idname, 'SPACE', 'PRESS', ctrl=True, shift=True)
    kmi.properties.total = 4
    addon_keymaps.append(km)

def unregister():
    bpy.utils.unregister_class(AddEdgeSlabs)
    bpy.types.VIEW3D_MT_object.remove(menu_func)

    # handle the keymap
    wm = bpy.context.window_manager
    for km in addon_keymaps:
        wm.keyconfigs.addon.keymaps.remove(km)
    # clear the list
    del addon_keymaps[:]


if __name__ == "__main__":
    register()
