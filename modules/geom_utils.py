"""Computational Geometry Utilities.

.. moduleauthor:: Nicola Cavallini

"""

import numpy as np
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import Polygon

from rtree import index

def map_tri_to_reference(tri,tri_map):
    tri = np.array(list(tri.exterior.coords))
    tri = tri[0:3,:].transpose()
    tri = np.vstack([tri,np.ones((1,3))])
    s_tri = np.dot(tri_map,tri)
    s_tri = s_tri[1:3,:].transpose()
    return s_tri


def tri_lin_map(xs_l,ys_l,s_l,t_l):
    A = np.vstack(
        [
        np.reshape(xs_l,(1,3)),
        np.reshape(ys_l,(1,3)),
        np.ones((1,3))
        ]
        )
    B = np.vstack(
        [
        np.reshape(s_l,(1,3)),
        np.reshape(t_l,(1,3)),
        np.ones((1,3))
        ]
        )

    tri_map = np.dot(B,np.linalg.inv(A))
    return tri_map

def triangulate(poly):
    # len(poly.exterior.coords)
    center = poly.centroid.coords[0]
    #print center
    n_edges = len(poly.exterior.coords)-1
    #print 'n_edges = ' + str(n_edges)
    p_list = poly.exterior.coords
    triangles = []
    for i in range(0,n_edges):
        p0 = p_list[i]
        p1 = p_list[i+1]
        #print p0
        tri = Polygon([center,p1,p0])
        triangles.append(tri)
        #print tri.area
    #print 
    #print 
    return triangles


def point_in_element(topo_u,x_u,y_u,XY_str,bbox_u):

    id_intersections = []

    for xy_s in XY_str:
        x_l = xy_s[0]
        y_l = xy_s[1]
        l = list(bbox_u.intersection( (x_l,y_l,x_l,y_l) ))
        id_intersections.append(l)

    point_in_tri = np.zeros(XY_str.shape[0], dtype=int)

    ipnt = 0
    for xy_s in XY_str:
        point = Point(tuple(xy_s.tolist()))
        tri_list = id_intersections[ipnt]
        for tri in tri_list:
            eval_p = np.zeros((3,2))
            eval_p[:,0] = x_u[topo_u[tri,:]]
            eval_p[:,1] = y_u[topo_u[tri,:]]
            el = Polygon(tuple(eval_p.tolist()))
            touch = el.touches(point)
            contains = el.contains(point)
            if touch == True or contains == True:
                point_in_tri[ipnt] = tri
        ipnt += 1
    
    return point_in_tri

def make_bounding_box(topo_u,x_u,y_u):
    bb_u = index.Index()
    iel = 0
    for row in topo_u:
        x_l = x_u[row]
        y_l = y_u[row]
        left, bottom, right, top = (min(x_l),min(y_l),max(x_l),max(y_l))
        bb_u.insert(iel,(left, bottom, right, top))
        iel += 1
    return bb_u

def fluid_intersect_string(topo_u,x_u,y_u,
                           topo_s,x_s,y_s):
    """
    Intersection between the fluid mesh and the solid mesh,
    just the same code a before but now the fluid is intersecting 
    a string.

    Input:
    
    topo_u,x_u,y_u : The fluid mesh.\n
    topo_s,x_s,y_s : The soolid mesh.\n
    
    Output:
    
    str_segments: each row of this array correspond to an original 
    solid mesh element. Each elememt on a row, is a chunk in wich
    the original solid mesh element
    is decomposed.\n
    
    fluid_id: the fluid element corresponding to each tructure chunk.
    
    """

    bb_u = make_bounding_box(topo_u,x_u,y_u)

    id_intersections = []

    for row in topo_s:
        x_l = x_s[row]
        y_l = y_s[row]
        l = list(bb_u.intersection((min(x_l),min(y_l),max(x_l),max(y_l))))
        id_intersections.append(l)

    #print id_intersections

    str_iel = 0
    str_segments = []
    fluid_id = []
    for row in id_intersections:
        x_l = x_s[topo_s[str_iel,:]]
        y_l = y_s[topo_s[str_iel,:]]
        eval_p = np.zeros((x_l.shape[0],2))
        eval_p[:,0] = x_l
        eval_p[:,1] = y_l
        line = LineString(tuple(eval_p.tolist()))
        #print line
        chunks = []
        chunk_el = []
        for el in row:
            x_l = x_u[topo_u[el,:]]
            y_l = y_u[topo_u[el,:]]
            eval_p = np.zeros((x_l.shape[0],2))
            eval_p[:,0] = x_l
            eval_p[:,1] = y_l
            polygon = Polygon(tuple(eval_p.tolist()))
            #query = line.intersects(polygon)
            intersection = line.intersection(polygon)
            #print polygon
            #print quer
            #print intersection.geom_type
            if intersection.geom_type == 'LineString':
                chunks.append(intersection)
                chunk_el.append(el)
                #print 'TROVATA!'
        str_segments.append(chunks)
        fluid_id.append(chunk_el)
        str_iel += 1
    return str_segments,fluid_id

def fluid_intersect_mesh(topo_u,x_u,y_u,
                         topo_s,x_s,y_s):
    """
    Intersection between the fluid mesh and the solid mesh.

    Input:
    
    topo_u,x_u,y_u : The fluid mesh.\n
    topo_s,x_s,y_s : The soolid mesh.\n
    
    Output:
    
    str_segments: each row of this array correspond to an original 
    solid mesh element. Each elememt on a row, is a chunk in wich
    the original solid mesh element
    is decomposed.\n
    
    fluid_id: the fluid element corresponding to each tructure chunk.
    
    """
    bb_u = make_bounding_box(topo_u,x_u,y_u)

    id_intersections = []

    for row in topo_s:
        x_l = x_s[row]
        y_l = y_s[row]
        l = list(bb_u.intersection((min(x_l),min(y_l),max(x_l),max(y_l))))
        id_intersections.append(l)

    str_iel = 0
    str_segments = []
    fluid_id = []
    for row in id_intersections:
        x_l = x_s[topo_s[str_iel,:]]
        y_l = y_s[topo_s[str_iel,:]]
        eval_p = np.zeros((x_l.shape[0],2))
        eval_p[:,0] = x_l
        eval_p[:,1] = y_l
        line = Polygon(tuple(eval_p.tolist()))
        #print line
        chunks = []
        chunk_el = []
        for el in row:
            x_l = x_u[topo_u[el,:]]
            y_l = y_u[topo_u[el,:]]
            eval_p = np.zeros((x_l.shape[0],2))
            eval_p[:,0] = x_l
            eval_p[:,1] = y_l
            polygon = Polygon(tuple(eval_p.tolist()))
            #query = line.intersects(polygon)
            intersection = line.intersection(polygon)
            #print polygon
            #print quer
            #print intersection.geom_type
            if intersection.geom_type == 'Polygon':
                chunks.append(intersection)
                chunk_el.append(el)
                #print 'TROVATA!'
        str_segments.append(chunks)
        fluid_id.append(chunk_el)
        str_iel += 1
    return str_segments,fluid_id

def get_reference_coords(topo_s,x_s,y_s,s_lgr,str_iel,coord_list):
    #
    # WARNING: 
    ds = s_lgr[1] - s_lgr[0]#THIS LINE ONLY WORKS FOR UNIFORM GRIDS
    #
    #
    s0 = s_lgr[topo_s[str_iel,0]]
    s1 = s_lgr[topo_s[str_iel,1]]
    x0 = x_s[topo_s[str_iel,0]]
    x1 = x_s[topo_s[str_iel,1]]
    y0 = y_s[topo_s[str_iel,0]]
    y1 = y_s[topo_s[str_iel,1]]
    se_l = LineString([(x0,y0),(x1,y1)])
    sp = np.zeros((len(coord_list),1))
    p_id = 0
    for point in coord_list:
        p = Point(point)
        s = s0 + np.sqrt((p.x-x0)**2+(p.y-y0)**2)/se_l.length * ds
        sp[p_id,0] = s 
        p_id += 1
        
    return sp

def area_evaluation(x,y):
    eval_p = np.zeros((x.shape[0],2))
    eval_p[:,0] = x
    eval_p[:,1] = y
    polygon = Polygon(tuple(eval_p.tolist()))
    return polygon.area