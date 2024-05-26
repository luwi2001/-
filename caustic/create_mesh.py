from utilities import Point3D,_centroid,Triangle,Mesh,_triangle_area
from PIL import Image
import math
import numpy as np


def squareMesh(width: int, height: int):
    nodeList = [Point3D() for _ in range(height * width)]
    nodeArray = [[Point3D() for _ in range(height)] for _ in range(width)]
    count = 0
    for y in range(height):
        for x in range(width):
            newPoint = Point3D(x, y, 0, x, y)
            nodeList[count] = newPoint
            nodeArray[x][y] = newPoint
            count += 1
    triangles = [Triangle() for _ in range((width - 1) * (height - 1) * 2)]
    count = 0
    for y in range(height-1):
        for x in range(width-1):
            # here x and y establish the column of squares we're in
            index_ul = y * width + x
            index_ur = index_ul + 1
            index_ll = (y + 1) * width + x
            index_lr = index_ll + 1
            triangles[count] = Triangle(index_ul, index_ll, index_ur)
            count += 1
            triangles[count] = Triangle(index_lr, index_ur, index_ll)
            count += 1
    return Mesh(nodeList, nodeArray, triangles, width, height)

def centroid(mesh: Mesh, index: int):
    triangle = mesh.triangles[index]
    p1 = mesh.nodes[triangle.pt1]
    p2 = mesh.nodes[triangle.pt2]
    p3 = mesh.nodes[triangle.pt3]
    return _centroid(p1, p2, p3)




def findT(p1, p2, p3, dp1, dp2, dp3):
    x1 = p2.x - p1.x
    y1 = p2.y - p1.y
    x2 = p3.x - p1.x
    y2 = p3.y - p1.y
    u1 = dp2.x - dp1.x
    v1 = dp2.y - dp1.y
    u2 = dp3.x - dp1.x
    v2 = dp3.y - dp1.y
    a = u1 * v2 - u2 * v1
    b = x1 * v1 + y2 * u1 - x2 * v1 - y1 * u2
    c = x1 * y2 - x2 * y1
    if a != 0:
        quotient = b**2 - 4*a * c
        if quotient >= 0:
            d = math.sqrt(quotient)
            return (-b - d) / (2*a), (-b + d) / (2*a)
        else:
            return -123.0, -123.0
    else:
        try:
            return -c / b, -c / b
        except(ZeroDivisionError):
            return -123.0, -123.0

def saveObj(mesh, filename, scale=1.0, scalez=1.0, reverse=False, flipxy=False):
    # This function saves the mesh object in stl format
    with open(filename, "w") as io:
        for vertex in mesh.nodes:
            if flipxy:
                io.write("v " + str(vertex.y * scale) + " " + str(vertex.x * scale) + " " + str(vertex.z * scalez) + "\n")
            else:
                io.write("v " + str(vertex.x * scale) + " " + str(vertex.y * scale) + " " + str(vertex.z * scalez) + "\n")
        for face in mesh.triangles:
            if reverse:
                io.write("f " + str(face.pt3) + " " + str(face.pt2) + " " + str(face.pt1) + "\n")
            else:
                io.write("f " + str(face.pt1) + " " + str(face.pt2) + " " + str(face.pt3) + "\n")
                
        
        io.write("dims " + str(mesh.width) + " " + str(mesh.height) + "\n")


def gradient(f):
    width, height = np.shape(f)
    delta_fu = np.zeros((width, height))
    delta_fv = np.zeros((width, height))
    for x in range(width):
        for y in range(height):
            delta_fu[x, y] = 0 if x == width-1 else f[x + 1, y] - f[x, y]
            delta_fv[x, y] = 0 if y == height-1 else f[x, y + 1] - f[x, y]
    return delta_fu, delta_fv



def getPixelArea(mesh: Mesh):
    pixelAreas = np.zeros((mesh.width-1, mesh.height-1))
    for x in range(mesh.width-1):
        for y in range(mesh.height-1):
            upperLeft = mesh.nodeArray[x][y]
            upperRight = mesh.nodeArray[x + 1][y]
            lowerLeft = mesh.nodeArray[x][y + 1]
            lowerRight = mesh.nodeArray[x + 1][y + 1]
            pixelAreas[x, y] = _triangle_area(lowerLeft, upperRight, upperLeft) + _triangle_area(lowerLeft, lowerRight, upperRight)
    return pixelAreas




def relax(grid):
    newgrid = grid.copy()
    
    for i in range(1,newgrid.shape[0]-1):
        for j in range(1,newgrid.shape[1]-1):
            newgrid[i,j] = 0.25 * (newgrid[i,j+1] + newgrid[i,j-1] +
                                   newgrid[i+1,j] + newgrid[i-1,j])
    
    return newgrid


def marchMesh(mesh, phi):
    grad_phi_u, grad_phi_v = np.gradient(phi)
    imgWidth, imgHeight = phi.shape
    velocities = np.empty((mesh.width, mesh.height), dtype=Point3D)
    for x in range(mesh.width):
        for y in range(mesh.height):
            if x == mesh.width - 1:
                u = 0
            else:
                u = grad_phi_u[x, y - 1] if y == mesh.height - 1 else grad_phi_u[x, y]
            if y == mesh.height - 1:
                v = 0
            else:
                v = grad_phi_v[x - 1, y] if x == mesh.width - 1 else grad_phi_v[x, y]
            velocities[x, y] = Point3D(-u, -v, 0, 0, 0)
    min_t = 10000
    for triangle in mesh.triangles:
            p1 = mesh.nodes[triangle.pt1]
            p2 = mesh.nodes[triangle.pt2]
            p3 = mesh.nodes[triangle.pt3]
            v1 = velocities[p1.ix, p1.iy]
            v2 = velocities[p2.ix, p2.iy]
            v3 = velocities[p3.ix, p3.iy]
            t1, t2 = findT(p1, p2, p3, v1, v2, v3)
            if 0 < t1 < min_t:
                min_t = t1
            if 0 < t2 < min_t:
                min_t = t2
    print("Overall min_t:", min_t) 
    delta = min_t / 2
    for point in mesh.nodes:
        v = velocities[point.ix, point.iy]
        point.x = v.x * delta + point.x
        point.y = v.y * delta + point.y
        
def quantifyLoss(D, suffix, img):
    print("Loss:")
    print("Minimum:", np.min(D))
    print("Maximum:", np.max(D))
    print("SUM:", np.sum(D))


def oneIteration(meshy, img, suffix):
    LJ = getPixelArea(meshy)
    D = np.float64(LJ - img)
    D = D - D.mean()
    print(np.min(D))
    print(np.max(D))
    quantifyLoss(D, suffix, img)
    width, height = img.shape
    ϕ = np.zeros((width, height))
    print("Building Phi")
    for i in range(1000):
        ϕ = relax(D)
    marchMesh(meshy, ϕ)

def setHeights(mesh, heights, heightScale=1.0, heightOffset=10):
    width, height = heights.shape
    for y in range(height):
        for x in range(width):
            mesh.nodeArray[x][y].z = heights[x, y] * heightScale + heightOffset
    for y in range(height+1):
        mesh.nodeArray[width][y].z = mesh.nodeArray[width-1][y].z
    for x in range(width+1):
        mesh.nodeArray[x][height].z = mesh.nodeArray[x][height-1].z

    print("We now have $(count - 1) valid nodes")



def solidify(inputMesh, offset=1):
    width = inputMesh.width
    height = inputMesh.height
    totalNodes = width * height * 2
    nodeList = np.empty(totalNodes, dtype=object)
    nodeArrayTop = np.empty((width, height), dtype=object)
    nodeArrayBottom = np.empty((width, height), dtype=object)
    numEdgeNodes = width * 2 + (height - 2) * 2
    numTrianglesTop = (width-1)  * (height-1) * 2
    numTrianglesBottom = numTrianglesTop
    numTrianglesEdges = numEdgeNodes * 2
    totalTriangles = numTrianglesBottom + numTrianglesTop + numTrianglesEdges
    print(f"Specs: {width}  {height}  {totalNodes}  {numEdgeNodes}  {numTrianglesBottom} {totalTriangles}")
    # Build the bottom surface
    count = 0
    for y in range(height):
        for x in range(width):
            newPoint = Point3D(x, y, -offset, x, y)
            nodeList[count] = newPoint
            nodeArrayBottom[x, y] = newPoint
            count += 1
            
    for y in range(height):
        for x in range(width):
            node = inputMesh.nodeArray[x][y]
            copiedPoint = Point3D(node.x, node.y, node.z, node.ix, node.iy)
            if node.ix != x:
               #print("OH NO POINTS x")
               print('ix:',node.ix,"        ","x:",x)
            if node.iy != y:
               print('iy:',node.iy,"        ","y:",y)
               #print("OH NO POINTS y")
            nodeList[count] = copiedPoint
            nodeArrayTop[x, y] = copiedPoint
            count += 1
    print(f"We now have {count} valid")

    triangles = [Triangle() for _ in range(totalTriangles)]
    count = 0
    
    
    
    for y in range(height-1):
        for x in range(width-1):
            index_ul = (y) * width + x
            index_ur = index_ul + 1
            index_ll = (y + 1) * width + x
            index_lr = index_ll + 1
            triangles[count] = Triangle(index_ul, index_ll, index_ur)
            count += 1
            triangles[count] = Triangle(index_lr, index_ur, index_ll)
            count += 1
    print(f"We've filled up {count} triangles")
    if count != numTrianglesBottom:
        print(f"Hmm aren't count and triangles bottom equal? {count} vs {numTrianglesBottom}")
    for y in range(height-1):
        for x in range(width-1):
            index_ul = y * width + x + totalNodes // 2
            index_ur = index_ul + 1
            index_ll = (y + 1) * width + x + totalNodes // 2
            index_lr = index_ll + 1
            triangles[count] = Triangle(index_ul, index_ur, index_ll)
            count += 1
            triangles[count] = Triangle(index_lr, index_ll, index_ur)
            count += 1
    print(f"We've filled up {count} triangles")
    
    
    ####edge triangles
    x = 0
    for y in range(height-1):
        ll = y * width + x
        ul = ll + totalNodes // 2
        lr = (y + 1) * width + x
        ur = lr + totalNodes // 2
        triangles[count] = Triangle(ll, ul, ur)
        count += 1
        triangles[count] = Triangle(ur, lr, ll)
        count += 1
    
    print(f"We've filled up {count} triangles")
    
    x = width - 1
    for y in range(height-1):
        ll = y * width + x
        ul = ll + totalNodes // 2
        lr = (y + 1) * width + x
        ur = lr + totalNodes // 2
        triangles[count] = Triangle(ll, ur, ul)
        count += 1
        triangles[count] = Triangle(ur, ll, lr)
        count += 1
    
    print(f"We've filled up {count} triangles")
    
    y = 0
    for x in range(1, width):
        ll = y * width + x
        ul = ll + totalNodes // 2
        lr = y * width + (x - 1)
        ur = lr + totalNodes // 2
        triangles[count] = Triangle(ll, ul, ur)
        count += 1
        triangles[count] = Triangle(ur, lr, ll)
        count += 1
    
    print(f"We've filled up {count} triangles")
    
    y = height - 1
    for x in range(1,width):
        ll = y * width + x
        ul = ll + totalNodes // 2
        lr = y * width + (x - 1)
        ur = lr + totalNodes // 2
        triangles[count] = Triangle(ll, ur, ul)
        count += 1
        triangles[count] = Triangle(ur, ll, lr)
        count += 1
    
    print(f"We've filled up {count} triangles")
    return Mesh(nodeList, nodeArrayBottom, triangles, width, height)

def solidify2(inputMesh, offset=1):
    width = inputMesh.width
    height = inputMesh.height
    totalNodes = width * height + 4
    nodeList = np.empty(totalNodes, dtype=object)
    nodeArrayTop = np.empty((width, height), dtype=object)
    nodeArrayBottom = np.empty((width, height), dtype=object)
    numTrianglesTop = (width-1)  * (height-1) * 2
    numTrianglesBottom = 2
    numTrianglesEdges = 4 * 2
    totalTriangles = numTrianglesBottom + numTrianglesTop + numTrianglesEdges
    #print(f"Specs: {width}  {height}  {totalNodes}  {numEdgeNodes}  {numTrianglesBottom} {totalTriangles}")
    # Build the bottom surface
    count = 0
    for y in [0, height-1]:
        for x in [0, width-1]:
            newPoint = Point3D(x, y, -offset, x, y)
            nodeList[count] = newPoint
            nodeArrayBottom[x, y] = newPoint
            count += 1
            
    for y in range(height):
        for x in range(width):
            node = inputMesh.nodeArray[x][y]
            copiedPoint = Point3D(node.x, node.y, node.z, node.ix, node.iy)
            if node.ix != x:
               #print("OH NO POINTS x")
               print('ix:',node.ix,"        ","x:",x)
            if node.iy != y:
               print('iy:',node.iy,"        ","y:",y)
               #print("OH NO POINTS y")
            nodeList[count] = copiedPoint
            nodeArrayTop[x, y] = copiedPoint
            count += 1
    print(f"We now have {count} valid")

    triangles = [Triangle() for _ in range(totalTriangles)]
    
    count = 0   
    
    index_ul = 0
    index_ur = 1
    index_ll = 2
    index_lr = 3
    triangles[count] = Triangle(index_ul, index_ll, index_ur)
    count += 1
    triangles[count] = Triangle(index_lr, index_ur, index_ll)
    count += 1
    
    print(f"We've filled up {count} triangles")
    if count != numTrianglesBottom:
        print(f"Hmm aren't count and triangles bottom equal? {count} vs {numTrianglesBottom}")
    for y in range(height-1):
        for x in range(width-1):
            index_ul = y * width + x + 4
            index_ur = index_ul + 1
            index_ll = (y + 1) * width + x + 4
            index_lr = index_ll + 1
            triangles[count] = Triangle(index_ul, index_ur, index_ll)
            count += 1
            triangles[count] = Triangle(index_lr, index_ll, index_ur)
            count += 1
    print(f"We've filled up {count} triangles")
    
    
    
    b1 = 0
    b2 = 1
    b3 = 2
    b4 = 3
    t1 = 4
    t2 = 4 + width - 1
    t3 = 4 + width * height - width
    t4 = 4 + width * height - 1
    
    triangles[count] = Triangle(b1, t1, t2)
    count += 1
    
    triangles[count] = Triangle(b1, b2, t2)
    count += 1    
    
    triangles[count] = Triangle(b2, t2, t4)
    count += 1    
    
    triangles[count] = Triangle(b2, b4, t4)
    count += 1    
    
    triangles[count] = Triangle(b4, t4, t3)
    count += 1    
    
    triangles[count] = Triangle(b4, b3, t3)
    count += 1   
 
    triangles[count] = Triangle(b3, t3, t1)
    count += 1
    
    triangles[count] = Triangle(b3, b1, t1)
    count += 1
    
    # ####edge triangles
    # ll = 0
    # ul = ll + 4
    # lr = 1
    # ur = 4 + width - 1
    # triangles[count] = Triangle(ll, ul, ur)
    # count += 1
    # triangles[count] = Triangle(ur, lr, ll)
    # count += 1
    
    # print(f"We've filled up {count} triangles")
               
    # ll = 1
    # ul = 4 + height - 1
    # lr = 3
    # ur =  width * height + 4 - 1
    # triangles[count] = Triangle(ll, ur, ul)
    # count += 1
    # triangles[count] = Triangle(ur, ll, lr)
    # count += 1
    
    # print(f"We've filled up {count} triangles")
    
        
    # ll = 3
    # ul = width * height + 4 - 1
    # lr = 2
    # ur = width * height + 4 - width 
    # triangles[count] = Triangle(ll, ul, ur)
    # count += 1
    # triangles[count] = Triangle(ur, lr, ll)
    # count += 1
    
    # print(f"We've filled up {count} triangles")
    
        
    # ll = 2
    # ul = width * height + 4 - width 
    # lr = 0
    # ur = lr + 4
    # triangles[count] = Triangle(ll, ur, ul)
    # count += 1
    # triangles[count] = Triangle(ur, ll, lr)
    # count += 1
    
    # print(f"We've filled up {count} triangles")
    return Mesh(nodeList, nodeArrayBottom, triangles, width, height)

def findSurface(mesh, image, f, imgWidth):
    width, height = image.shape
    H = f
    metersPerPixel = imgWidth / width
    print(metersPerPixel)
    n2 = 1
    n1 = 1.49
    Nx = np.zeros((width + 1, height + 1))
    Ny = np.zeros((width + 1, height + 1))
    for j in range(height):
        for i in range(width):
            node = mesh.nodeArray[i][j]
            dx = (node.ix - node.x) * metersPerPixel
            dy = (node.iy - node.y) * metersPerPixel
            little_h = node.z * metersPerPixel
            H_minus_h = H - little_h
            dz = H_minus_h
            Ny[i, j] = np.tan(np.arctan(dy / dz) / (n1 - n2))
            Nx[i, j] = np.tan(np.arctan(dx / dz) / (n1 - n2))
    divergence = np.zeros((width, height))
    for j in range(height):
        for i in range(width):
            δx = (Nx[i + 1, j] - Nx[i, j])
            δy = (Ny[i, j + 1] - Ny[i, j])
            divergence[i, j] = δx + δy
    print("Have all the divergences")
    print("Divergence sum: {}".format(np.sum(divergence)))
    divergence = divergence - divergence.mean()
    
    h = np.zeros((width, height))
    for i in range(2000):
        h = relax(divergence)
    return h, metersPerPixel




def engineer_caustics(img):
    img = np.array(img, dtype=np.uint8)
    img2 = np.transpose(img) * 1
    width, height = img2.shape
    meshy = squareMesh(width + 1, height + 1)
    mesh_sum = width * height
    image_sum = np.sum(img2)
    boost_ratio = mesh_sum / image_sum
    img3 = img2 * boost_ratio
    for i in range(3):
        oneIteration(meshy, img3, "it1")
    artifactSize = 0.1  # meters
    focalLength = 0.2 # meters
    h, metersPerPixel = findSurface(meshy, img3, focalLength, artifactSize)
    setHeights(meshy, h, 1.0, 0.01)
    solidMesh = solidify(meshy)
    return meshy,solidMesh,img3




# def main():
#     img = Image.open('dog_2.png')
#     img = img.convert('L')
#     return engineer_caustics(img)

if __name__ == '__main__':
    #meshy,solidMesh, img3 = main()
    img = Image.open('dog_2.png')
    img = img.convert('L')
    img = np.array(img, dtype=np.uint8)
    img2 = np.transpose(img) * 1
    width, height = img2.shape
    meshy = squareMesh(width + 1, height + 1)
    mesh_sum = width * height
    image_sum = np.sum(img2)
    boost_ratio = mesh_sum / image_sum
    img3 = img2 * boost_ratio
    for i in range(3):
        oneIteration(meshy, img3, "it1")
    artifactSize = 0.1  # meters
    focalLength = 0.2 # meters
    h, metersPerPixel = findSurface(meshy, img3, focalLength, artifactSize)
    meshy_1 = meshy
    setHeights(meshy_1, h, 1.0, 10)
    solidMesh = solidify2(meshy_1)
    
    # 
    import open3d as o3d
    
    vertices= []
    for node in solidMesh.nodes:
        node_coor = [node.x/solidMesh.width,node.y/solidMesh.height,max(node.z/solidMesh.width*300-11.3,0)]
        vertices.append(node_coor)
    vertices = np.array(vertices)
    vertices = vertices.round(4)
    
  
    face=[]
    for tri in solidMesh.triangles:
        face_coor = (tri.pt1,tri.pt2,tri.pt3)
        face.append(face_coor)
    face = np.array(face)

    def get_non_manifold_vertex_mesh(verts, triangles):
    
        mesh = o3d.geometry.TriangleMesh()
        mesh.vertices = o3d.utility.Vector3dVector(verts)
        mesh.triangles = o3d.utility.Vector3iVector(triangles)
        #mesh.compute_vertex_normals()  # 估计法向量后支持光照渲染
        
        mesh.compute_triangle_normals()

        
        return mesh

    mesh_out = get_non_manifold_vertex_mesh(vertices,face)
  
    
    o3d.io.write_triangle_mesh(filename="obj_file/out2.obj",
                               mesh=mesh_out,
                               write_ascii=False,
                               compressed=False,
                               write_vertex_normals=False,
                               write_vertex_colors=False,
                               write_triangle_uvs=False,
                               print_progress=True,
                               )
    
    
