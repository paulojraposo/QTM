# -*- coding: utf-8 -*-

#   .-.
#   /v\    L   I   N   U   X
#  // \\
# /(   )\
#  ^^-^^

# This script makes a Quarternary Triangular Mesh (QTM) to tessellate the world based
# on an octohedron, based on Geoffrey Dutton's conception:
#
# Dutton, Geoffrey H. "Planetary Modelling via Hierarchical Tessellation." In Procedings of the
# AutoCarto 9 Conference, 462–71. Baltimore, MD, 1989.
# https://pdfs.semanticscholar.org/875e/12ce948012b3eced58c0f1470dba51ef87ef.pdf
#
# This script written by Paulo Raposo (pauloj.raposo [at] outlook.com) and Randall Brown
# (ranbrown8448 [at] gmail.com). Under MIT license.
#
# Dependencies: nvector (see https://pypi.python.org/pypi/nvector and http://www.navlab.net/nvector),
#               OGR Python bindings (packaged with GDAL).
#


import os, argparse, logging, datetime, sys
import nvector as nv
from osgeo import ogr, osr


def GetGeodeticMidpoint(vert1, vert2):
    """Given two Vertices, return the geodetic midpoint of the great circle arc between them, on the WGS84 ellipsoid. Uses nvector."""
    # see http://nvector.readthedocs.org/en/latest/src/overview.html?highlight=midpoint#description
    wgs84 = nv.FrameE(name='WGS84')
    n_EB_E_t0 = wgs84.GeoPoint(vert1[0], vert1[1], degrees=True).to_nvector()
    n_EB_E_t1 = wgs84.GeoPoint(vert2[0], vert2[1], degrees=True).to_nvector()
    path = nv.GeoPath(n_EB_E_t0, n_EB_E_t1)
    halfway = 0.5
    g_EB_E_ti = path.interpolate(halfway).to_geo_point()
    lat_ti, lon_ti = g_EB_E_ti.latitude_deg, g_EB_E_ti.longitude_deg
    return (float(lat_ti), float(lon_ti))


def constructGeometry(facet):
    """Accepting a list from this script that stores vertices, return an OGR Geometry polygon object."""
    ring = ogr.Geometry(ogr.wkbLinearRing)
    if len(facet) == 5:
        # This is a triangle facet of format (vert,vert,vert,vert,orient)
        vertexTuples = facet[:4]
    if len(facet) == 6:
        # This is a rectangle facet of format (vert,vert,vert,vert,vert,northboolean)
        vertexTuples = facet[:5]
    for vT in vertexTuples:
        ring.AddPoint(vT[1], vT[0]) # sequence: lon, lat (x,y)
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


def divideFacet(aFacet):
    """Will always return four facets, given one, rectangle or triangle."""

    # Important: For all facets, first vertex built is always the most south-then-west, going counter-clockwise thereafter.

    if len(aFacet) == 5:

        # This is a triangle facet.

        orient = aFacet[4]  # get the string expressing this triangle's orientation

        #       Cases, each needing subdivision:
        #                    ______           ___   ___
        #       |\      /|   \    /   /\     |  /   \  |     ^
        #       | \    / |    \  /   /  \    | /     \ |     N
        #       |__\  /__|     \/   /____\   |/       \|
        #
        #        up    up     down    up     down    down   -- orientations, as "u" or "d" in code below.


        # Find the geodetic bisectors of the three sides, store in sequence using edges defined
        # by aFacet vertex indeces: [0]&[1] , [1]&[2] , [2]&[3]
        newVerts = []
        for i in range(3):
            newVerts.append(GetGeodeticMidpoint(aFacet[i], aFacet[i + 1]))

        if orient == "u":
            #          In the case of up facets, there will be one "top" facet
            #          and 3 "bottom" facets after subdivision; we build them in the sequence inside the triangles:
            #
            #                   2
            #                  /\         Outside the triangle, a number is the index of the vertex in aFacet,
            #                 / 1\        and a number with an asterisk is the index of the vertex in newVerts.
            #             2* /____\ 1*
            #               /\ 0  /\
            #              /2 \  /3 \
            #             /____\/____\
            #           0or3   0*     1

            newFacet0 = [newVerts[0], newVerts[1], newVerts[2], newVerts[0], "d"]
            newFacet1 = [newVerts[2], newVerts[1], aFacet[2], newVerts[2], "u"]
            newFacet2 = [aFacet[0], newVerts[0], newVerts[2], aFacet[0], "u"]
            newFacet3 = [newVerts[0], aFacet[1], newVerts[1], newVerts[0], "u"]

        if orient == "d":
            #          In the case of down facets, there will be three "top" facets
            #          and 1 "bottom" facet after subdivision; we build them in the sequence inside the triangles:
            #
            #            2_____1*_____1
            #             \ 2  /\ 3  /
            #              \  / 0\  /    Outside the triangle, a number is the index of the vertex in aFacet,
            #               \/____\/     and a number with an asterisk is the index of the vertex in newVerts.
            #              2*\ 1  /0*
            #                 \  /
            #                  \/
            #                 0or3

            newFacet0 = [newVerts[2], newVerts[0], newVerts[1], newVerts[2], "u"]
            newFacet1 = [aFacet[0], newVerts[0], newVerts[2], aFacet[0], "d"]
            newFacet2 = [newVerts[2], newVerts[1], aFacet[2], newVerts[2], "d"]
            newFacet3 = [newVerts[0], aFacet[1], newVerts[1], newVerts[0], "d"]

    if len(aFacet) == 6:

        # This is a rectangle facet.

        northBoolean = aFacet[5]  # true for north, false for south

        if northBoolean:

            # North pole rectangular facet.

            # Build new facets in the sequence inside the polygons:

            #          3..........2   <-- North Pole
            #           |        |
            #           |   1    |    Outside the polys, a number is the index of the vertex in aFacet,
            #           |        |    and a number with an asterisk is the index of the vertex in newVerts.
            #           |        |
            #         2*|--------|1*           /\
            #           |\      /|  on globe  /__\
            #           | \ 0  / |  -------> /\  /\
            #           |  \  /  |          /__\/__\
            #           | 2 \/ 3 |
            #       0or4''''''''''1
            #               0*

            newVerts = []

            for i in range(4):
                if i != 2:
                    # on iter == 1 we're going across the north pole - don't need this midpoint.
                    newVerts.append(GetGeodeticMidpoint(aFacet[i], aFacet[i + 1]))

            newFacet0 = [newVerts[0], newVerts[1], newVerts[2], newVerts[0], "d"]  # triangle
            newFacet1 = [newVerts[2], newVerts[1], aFacet[2], aFacet[3], newVerts[2], True]  # rectangle
            newFacet2 = [aFacet[0], newVerts[0], newVerts[2], aFacet[0], "u"]  # triangle
            newFacet3 = [newVerts[0], aFacet[1], newVerts[1], newVerts[0], "u"]  # triangle

        else:

            # South pole rectangular facet

            #               1*
            #          3..........2
            #           | 2 /\ 3 |     Outside the polys, a number is the index of the vertex in aFacet,
            #           |  /  \  |     and a number with an asterisk is the index of the vertex in newVerts.
            #           | / 0  \ |
            #           |/      \|           ________
            #         2*|--------|0*         \  /\  /
            #           |        |  on globe  \/__\/
            #           |   1    |  ------->   \  /
            #           |        |              \/
            #           |        |
            #       0or4'''''''''1   <-- South Pole

            newVerts = []

            for i in range(4):
                if i != 0:
                    # on iter == 3 we're going across the south pole - don't need this midpoint
                    newVerts.append(GetGeodeticMidpoint(aFacet[i], aFacet[i + 1]))

            newFacet0 = [newVerts[2], newVerts[0], newVerts[1], newVerts[2], "u"]  # triangle
            newFacet1 = [aFacet[0], aFacet[1], newVerts[0], newVerts[2], aFacet[0], False]  # rectangle
            newFacet2 = [newVerts[2], newVerts[1], aFacet[3], newVerts[2], "d"]  # triangle
            newFacet3 = [newVerts[1], newVerts[0], aFacet[2], newVerts[1], "d"]  # triangle


    # In all cases, return the four facets made in a list
    return [newFacet0, newFacet1, newFacet2, newFacet3]


def printandlog(msg):
    """Given a string, this will both log it and print it to the console."""
    print(msg)
    logging.info(msg)


def main():
    # Input shell arguments
    parser = argparse.ArgumentParser(description='Builds a Dutton QTM (see citations in source code) and outputs it as a GeoJSON file in WGS84 coordinates.')

    # parser.add_argument('LEVELS', help='Integer number of levels to subdivide the QTM. Minimum of 1.')
    parser.add_argument('OUTSHPFILEDIR', help='Full path to output directory for the product QTM shapefiles.')
    parser.add_argument('LEVELS', help='Number of levels to generate. Give as an integer.')
    args = parser.parse_args()

    nLevels = int(args.LEVELS)
    outFileDir = args.OUTSHPFILEDIR

    # Log file setup
    dirPath = outFileDir
    logFile = os.path.join(dirPath, "qtm_creation_log.txt")
    logging.basicConfig(filename=logFile, level=logging.DEBUG)

    startTime = datetime.datetime.now()
    printandlog("Starting, " + str(startTime))

    printandlog("Total levels requested: " + str(nLevels))

    # The following WKT string from http://spatialreference.org/ref/epsg/4326/
    wktCoordSys = """GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]"""

    # Build the vertices of the initial 8 octohedron facets, as rectangles. Tuples: (Lat,Lon)
    #
    # N Pole  ---> (90,-180) ------- (90,-90)  ------- (90,0)  -------- (90, 90)  ------ (90, 180)
    #                  |                 |                |                |                 |
    #                  |                 |          Prime Meridian         |                 |
    #                  |                 |                |                |                 |
    # Equator ---> (0, -180)  ------ (0, -90)  ------- (0,0)  ---------  (0, 90)  ------ (0, 180)
    #                  |                 |                |                |                 |
    #                  |                 |                |                |                 |
    #                  |                 |                |                |                 |
    # S Pole  ---> (-90, -180) ----- (-90, -90) ------ (-90,0) -------- (-90, 90) ------ (-90, 180)

    p90_n180 = (90.0, -180.0)
    p90_n90 = (90.0, -90.0)
    p90_p0 = (90.0, 0.0)
    p90_p90 = (90.0, 90.0)
    p90_p180 = (90.0, 180.0)
    p0_n180 = (0.0, -180.0)
    p0_n90 = (0.0, -90.0)
    p0_p0 = (0.0, 0.0)
    p0_p90 = (0.0, 90.0)
    p0_p180 = (0.0, 180.0)
    n90_n180 = (-90.0, -180.0)
    n90_n90 = (-90.0, -90.0)
    n90_p0 = (-90.0, 0.0)
    n90_p90 = (-90.0, 90.0)
    n90_p180 = (-90.0, 180.0)

    # Keeping track of levels
    levelShapefiles = []
    levelFacets = {}
    ID = {}
    for lvl in range(nLevels):

        print("")
        printandlog("Working on level " + str(lvl))

        levelFacets[lvl] = []
        ID[lvl] = []
        previousLevel = None

        # Prepare a shapefile for this level
        outSHPFileName = "qtmlvl" + str(lvl) + ".geojson"
        levelShapefiles.append(outSHPFileName)
        sRef = osr.SpatialReference()
        sRef.ImportFromWkt(wktCoordSys)
        driver = ogr.GetDriverByName('GeoJSON')
        outFile = os.path.join(outFileDir, outSHPFileName)
        dst_ds = driver.CreateDataSource(outFile)
        fName = os.path.splitext(os.path.split(outFile)[1])[0]
        dst_layer = dst_ds.CreateLayer(fName, sRef, geom_type=ogr.wkbPolygon)
        levelFieldName = 'ID'
        layer_defn = dst_layer.GetLayerDefn()
        new_field = ogr.FieldDefn(levelFieldName, ogr.OFTInteger)
        dst_layer.CreateField(new_field)

        if lvl == 0:

            # Need to build the first level from scratch - all rectangle facets.
            # Important: For all facets, first vertex is always the most south-then-west, going counter-clockwise thereafter.

            # northern hemisphere
            levelFacets[0].append([p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True])
            levelFacets[0].append([p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True])
            levelFacets[0].append([p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True])
            levelFacets[0].append([p0_p90, p0_p180,  p90_p180, p90_p90, p0_p90, True])
            # southern hemisphere
            levelFacets[0].append([n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False])
            levelFacets[0].append([n90_n90, n90_p0,  p0_p0, p0_n90, n90_n90, False])
            levelFacets[0].append([n90_p0, n90_p90, p0_p90,  p0_p0, n90_p0, False])
            levelFacets[0].append([n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False])

            i = 0
            for f in levelFacets[lvl]:
                ID[0].append(str(i + 1))
                feature = ogr.Feature(layer_defn)
                feature.SetField('ID', ID[0][i])
                facetGeometry = constructGeometry(f)
                feature.SetGeometry(facetGeometry)
                dst_layer.CreateFeature(feature)
                feature.Destroy()  # Destroy the feature to free resources
                i = i + 1

        else:

            # Build further levels by subdividing the facets from the previous level.

            previousLevel = lvl - 1
            previousFacets = levelFacets[previousLevel]
            previousId = ID[previousLevel]

            nToSubdivide = len(previousFacets)
            iterlabel = 1
            i = 0
            k = 0
            for pf in previousFacets:

                sys.stdout.flush()  # for progress messages on console

                theseFacets = divideFacet(pf)
                j = 0

                for tF in theseFacets:
                    # Write to this level's shapefile
                    ID[lvl].append(previousId[i] + str(j))
                    feature = ogr.Feature(layer_defn)
                    feature.SetField('ID', ID[lvl][k])
                    facetGeometry = constructGeometry(tF)
                    feature.SetGeometry(facetGeometry)
                    dst_layer.CreateFeature(feature)
                    feature.Destroy()  # Destroy the feature to free resources
                    # Keep track of facet info
                    levelFacets[lvl].append(tF)
                    j = j + 1
                    k = k + 1

                # for progress messages on console
                prcnt = round((float(iterlabel) / float(nToSubdivide)) * 100, 3)
                sys.stdout.write("\r")
                sys.stdout.write("Progress: " + str(iterlabel) + " of " + str(nToSubdivide) + " | " + str(prcnt) + r" %... | Elapsed: " + str(datetime.datetime.now() - startTime))
                iterlabel += 1
                i = i + 1

        if previousLevel:
            del levelFacets[previousLevel]  # to free resources

        dst_ds.Destroy()  # Destroy the data source to free resouces

    endTime = datetime.datetime.now()
    print("")
    printandlog("Finished, " + str(endTime))
    elapsed = endTime - startTime
    printandlog("Total time for " + str(nLevels) + " levels: " + str(elapsed))

    # fin
    exit()


if __name__ == '__main__':
    main()
