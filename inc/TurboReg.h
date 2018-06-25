typedef enum Transformation {
/*********************************************************************
 Generic geometric transformation.
 ********************************************************************/
GENERIC_TRANSFORMATION = -1,

/*********************************************************************
 A translation is described by a single point. It keeps area, angle,
 and orientation. A translation is determined by 2 parameters.
 ********************************************************************/
TRANSLATION = 2,

/*********************************************************************
 A single points determines the translation component of a rigid-body
 transformation. As the rotation is given by a scalar number, it is
 not natural to represent this number by coordinates of a point. The
 rigid-body transformation is determined by 3 parameters.
 ********************************************************************/
RIGID_BODY = 3,

/*********************************************************************
 A pair of points determines the combination of a translation, of
 a rotation, and of an isotropic scaling. Angles are conserved. A
 scaled rotation is determined by 4 parameters.
 ********************************************************************/
SCALED_ROTATION = 4,

/*********************************************************************
 Three points generate an affine transformation, which is any
 combination of translation, rotation, isotropic scaling, anisotropic
 scaling, shearing, and skewing. An affine transformation maps
 parallel lines onto parallel lines and is determined by 6 parameters.
 ********************************************************************/
AFFINE = 6,


/*********************************************************************
 Four points describe a bilinear transformation, where a point of
 coordinates (x, y) is mapped on a point of coordinates (u, v) such
 that u = p0 + p1 x + p2 y + p3 x y and v = q0 + q1 x + q2 y + q3 x y.
 Thus, u and v are both linear in x, and in y as well. The bilinear
 transformation is determined by 8 parameters.
 ********************************************************************/
BILINEAR = 8,

/*********************************************************************
 Minimal linear dimension of an image in the multiresolution pyramid.
 ********************************************************************/
MIN_SIZE = 12,



} Transformation;