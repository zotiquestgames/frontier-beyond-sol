#!/usr/bin/env python
"""
GMLMaker converts star data to yEd GML format.

Input is a *.csv file with the following fields:
HABHYG: HabHYG catalog number
HIP: Hipparcos catalog number
IS_HAB: is star habitable? (i.e., in HabCat database?)
DISPLAY_NAME: showy name to use for human readable map/list/whatever
HYG: HYG catalog number
BAYERFLAMSTEED: Bayer-Flamsteed name
GLIESE: gliese catalog number
BD
HD,
HR,
PROPER_NAME,
SPECTRAL_CLASS,
DISTANCE,
XG,
YG,
ZG,
ABSMAG

Requires Python 2.3 or later (uses CSV and optparse module)
"""

#******************************************************************************
# Imports
import csv
import math
import optparse

#******************************************************************************
# Class definitions

class LinkAlgorithm:
    """used to select line-determining algorithm"""
    CLOSEST_NEIGHBOR,\
    TWO_CLOSEST_NEIGHBORS,\
    STAR_INTERFERENCE,\
    MAX_AND_MIN,\
    MATCH_LUMINOSITY,\
    DIST_LIMIT,\
    DIST_LIMIT_AND_NEIGHBOR,\
    DIST_LIMIT_AND_TWO_NEIGHBOR,\
    CLOSEST_HABITABLE,\
    TWO_CLOSEST_HABITABLE,\
    DIST_LIMIT_AND_HABITABLE,\
    DIST_LIMIT_AND_TWO_HABITABLE,\
    TWO_CLOSEST_HABITABLE_AND_TWO_NEIGHBOR\
    = range(13)
    DEFAULT_ALGORITHM = TWO_CLOSEST_NEIGHBORS
    
class SpectralClass:
    """used to encode spectral class"""
    CLASS_O,\
    CLASS_B,\
    CLASS_A,\
    CLASS_F,\
    CLASS_G,\
    CLASS_K,\
    CLASS_M,\
    CLASS_X,\
    NUM_CLASSES\
    = range(9)
    
    encode = {
    CLASS_O:"O",\
    CLASS_B:"B",\
    CLASS_A:"A",\
    CLASS_F:"F",\
    CLASS_G:"G",\
    CLASS_K:"K",\
    CLASS_M:"M",\
    CLASS_X:"X"}
    
    decode = {\
    "O":CLASS_O,\
    "B":CLASS_B,\
    "A":CLASS_A,\
    "F":CLASS_F,\
    "G":CLASS_G,\
    "K":CLASS_K,\
    "M":CLASS_M,\
    "X":CLASS_X}
    
class StarNodeColors:
    """color used to display spectral class in *.gml file"""
    color = {
    SpectralClass.CLASS_O:"#FF99FF",\
    SpectralClass.CLASS_B:"#99FFFF",\
    SpectralClass.CLASS_A:"#00FFFF",\
    SpectralClass.CLASS_F:"#00FF00",\
    SpectralClass.CLASS_G:"#FFFF99",\
    SpectralClass.CLASS_K:"#FF9F40",\
    SpectralClass.CLASS_M:"#FF9999",\
    SpectralClass.CLASS_X:"#FFBBBB"}
    
    HABSTAR_COLOR = "#FFFFFF"
    HABSTAR_BORDER_COLOR = "#FF00FF"
    
class InputFields:
    """used to index into fields of input CSV file"""
    HABHYG,\
    HIP,\
    IS_HAB,\
    DISPLAY_NAME,\
    HYG,\
    BAYERFLAMSTEED,\
    GLIESE,\
    BD,\
    HD,\
    HR,\
    PROPER_NAME,\
    SPECTRAL_CLASS,\
    DISTANCE,\
    XG,\
    YG,\
    ZG,\
    ABS_MAG\
    = range(17)
    
class StarField:
    """used to index into a starNode tuple"""
    HABHYG,\
    DISPLAY_NAME,\
    SPECTRAL_CODE,\
    IS_HAB,\
    LUMINOSITY,\
    DISTANCE,\
    LOCALE\
    = range(7)

class Coord:
    """used to index into an XYZ locale tuple"""
    XG,\
    YG,\
    ZG\
    = range(3)

class BBCoord:
    """used to index into an bounding box"""
    LEFT_X,\
    TOP_Y,\
    RIGHT_X,\
    BOTTOM_Y\
    = range(4)
    
class SystemField:
    """used to index into a starSystem node"""
    INDEX,\
    SYSTEM_NAME,\
    BRIGHTEST_SPECTRALCLASS,\
    BRIGHTEST_NAME,\
    LOCALE,\
    DISTANCE,\
    SYSTEM_STARS,\
    AT_LEAST_ONE_HABITABLE,\
    PRINT_LOCALE\
    = range(9)
    
class LinkField:
    """used to index into a link node"""
    START_INDEX,\
    DEST_INDEX,\
    DISTANCE\
    = range(3)
   
#******************************************************************************
# Function definitions

#--------------------------------------    
def GetSystemAt(starSystems, locale):
    """return the index of any exisiting system within 0.3 parsecs
    (one light year) of locale. Stars closer than that are assumed
    to be part of a multiple star system."""
    for systemIndex in starSystems:
        system = starSystems[systemIndex]
        systemLocale = system[SystemField.LOCALE]
        if math.fabs(systemLocale[Coord.XG] - locale[Coord.XG]) < 0.3 and\
           math.fabs(systemLocale[Coord.YG] - locale[Coord.YG]) < 0.3 and\
           math.fabs(systemLocale[Coord.ZG] - locale[Coord.ZG]) < 0.3:
            return system[SystemField.INDEX]
    return None

#--------------------------------------    
def LoadStars(theFileName):
    """read star data *.csv file and load data into starSystems"""    
    starSystems = {}
    boundingBox = [9999.0, 9999.0, -9999.0, -9999.0]
    index = -1
    starSystemIndex = -1
    excludeHeaderLine = True
    try:
        reader = csv.reader(file(theFileName + ".csv"))
        for row in reader:
            if excludeHeaderLine:
                excludeHeaderLine = False
                continue
            index += 1
            # create starNode from row
            habHYG = int(row[InputFields.HABHYG])
            displayName = row[InputFields.DISPLAY_NAME]
            theSpectralClass = row[InputFields.SPECTRAL_CLASS]
            theSpectralClass = theSpectralClass[0:1]
            try:
                spectralCode = SpectralClass.decode[theSpectralClass]
            except KeyError:
                spectralCode = SpectralClass.CLASS_X
            if row[InputFields.IS_HAB]:
                isHab = True
            else:
                isHab = False
            if float(row[InputFields.ABS_MAG]) > 0:
                luminosity = float(row[InputFields.ABS_MAG])**3.5
            else:
                luminosity = 0.0
            distance = float(row[InputFields.DISTANCE])
            locale = (float(row[InputFields.XG]),\
                       float(row[InputFields.YG]),\
                       float(row[InputFields.ZG]))
            
            starNode = (habHYG, displayName, spectralCode, isHab,\
                        luminosity, distance, locale)
                        
            if locale[Coord.XG] < boundingBox[BBCoord.LEFT_X]:
                boundingBox[BBCoord.LEFT_X] = locale[Coord.XG]
            if locale[Coord.XG] > boundingBox[BBCoord.RIGHT_X]:
                boundingBox[BBCoord.RIGHT_X] = locale[Coord.XG]
            if locale[Coord.YG] < boundingBox[BBCoord.TOP_Y]:
                boundingBox[BBCoord.TOP_Y] = locale[Coord.YG]
            if locale[Coord.YG] > boundingBox[BBCoord.BOTTOM_Y]:
                boundingBox[BBCoord.BOTTOM_Y] = locale[Coord.YG]
            
            print "Loading #%d: %s" % (index, starNode[1])
            
            # find star system that the new star belongs in
            theSystemIndex = GetSystemAt(starSystems, locale)
            if theSystemIndex == None:  # create one
                starSystemIndex += 1
                theSystemIndex = starSystemIndex
                brightestSpectralClassCode = SpectralClass.CLASS_X
                systemStars = []
                starSystemNode = [starSystemIndex, "", \
                    brightestSpectralClassCode, displayName, locale,\
                    distance, systemStars, False, locale]
                starSystems[starSystemIndex] = starSystemNode # add it
                
            
            # Get star system    
            starSystemNode = starSystems[theSystemIndex]
            
            # Add the new star to the star system
            if isHab:
                starSystemNode[SystemField.AT_LEAST_ONE_HABITABLE] = True
            if len(starSystemNode[SystemField.SYSTEM_NAME]) > 0:
                starSystemNode[SystemField.SYSTEM_NAME] += ", "
            starSystemNode[SystemField.SYSTEM_NAME] += displayName
            if starSystemNode[SystemField.BRIGHTEST_SPECTRALCLASS] > spectralCode:
                starSystemNode[SystemField.BRIGHTEST_SPECTRALCLASS] = spectralCode
                starSystemNode[SystemField.BRIGHTEST_NAME] = displayName
            starSystemNode[SystemField.SYSTEM_STARS].append(starNode)   
                
            # store changes
            starSystems[theSystemIndex] =  starSystemNode   
            
    except csv.Error:
        print "LoadStars: Error parsing" + theFileName + ".csv"
        
    boundingBox[BBCoord.LEFT_X] = math.floor(boundingBox[BBCoord.LEFT_X])    
    boundingBox[BBCoord.TOP_Y] = math.floor(boundingBox[BBCoord.TOP_Y])    
    boundingBox[BBCoord.RIGHT_X] = math.ceil(boundingBox[BBCoord.RIGHT_X])    
    boundingBox[BBCoord.BOTTOM_Y] = math.ceil(boundingBox[BBCoord.BOTTOM_Y])    
    
    return starSystems, boundingBox
        
#--------------------------------------    
def TrueDistance(pointA, pointB):
    """given two cartesian co-ords, calculatate distance between them"""
    deltaX = pointA[Coord.XG] - pointB[Coord.XG]
    deltaY = pointA[Coord.YG] - pointB[Coord.YG]
    deltaZ = pointA[Coord.ZG] - pointB[Coord.ZG]
    return math.sqrt((deltaX * deltaX) + 
                        (deltaY * deltaY) + 
                        (deltaZ * deltaZ))

#--------------------------------------    
def AddLink(starLinks, linkFile, startIndex, startName, 
            destIndex, destName, linkLength):
    """insert new link into starLinks and echo to linkFile"""
    # Simulate C++ trinary operator x ? y : z
    # with Python [z, y][bool(x)]
    # theKey = smallerIndex:largerIndex
    # so if link x to y is stored, new link y to x will NOT be stored
    theKey = "%d:%d" % ([destIndex, startIndex][startIndex < destIndex],
                        [destIndex, startIndex][startIndex >= destIndex])
    if not starLinks.has_key(theKey):
        starLinks[theKey] = (startIndex, destIndex, linkLength)
        linkFile.write("[%s]\t%.1fpc to\t[%s]\n" % (startName, linkLength, destName))
    return starLinks
    
#--------------------------------------    
def LinkNeighbor_(aStarSystem, starSystems, starLinks, linkFile, 
                    forbiddenStars, onlyHab):
    """Create link between aStarSystem and closest qualified
    starSystem and insert into starLinks.  forbiddenStars is a list of 
    unqualified systems, and onlyHab flag commands that systems on both ends 
    of the link are habitable.
    Selected system is inserted into forbiddenStars so that this
    function can be called multiple times.
    
    aStarSystem = origin star
    starSystems = list of all stars
    starLinks = list of star links
    linkFile = link textfile
    forbiddenStars = list of stars not to link to
    onlyHab = "only link to habitable stars" flag"""
    if aStarSystem[SystemField.INDEX] not in forbiddenStars:
        forbiddenStars.append(aStarSystem[SystemField.INDEX])
    # if onlyHab specifies, omit non-hab origin stars
    if (not onlyHab) or (aStarSystem[SystemField.AT_LEAST_ONE_HABITABLE]):
        closestNeighborDistance = 9999.0
        closestNeighborIndex = None
        for systemIndex in starSystems:
            # do not link a forbidden system system
            # (includes origin system and systems selected in
            # prior passes through this function)
            if systemIndex not in forbiddenStars:
                destSystem = starSystems[systemIndex]
                # if onlyHab specifies, omit non-hab destination stars
                if (not onlyHab) or \
                    (destSystem[SystemField.AT_LEAST_ONE_HABITABLE]):
                    distance = TrueDistance(aStarSystem[SystemField.LOCALE], 
                                            destSystem[SystemField.LOCALE])
                    if distance < closestNeighborDistance:
                        closestNeighborDistance = distance
                        closestNeighborIndex = systemIndex
        if closestNeighborIndex != None:
            destSystem = starSystems[closestNeighborIndex]
            forbiddenStars.append(closestNeighborIndex)
            starLinks = AddLink(starLinks, linkFile, 
                    aStarSystem[SystemField.INDEX],
                    aStarSystem[SystemField.SYSTEM_NAME],
                    destSystem[SystemField.INDEX],
                    destSystem[SystemField.SYSTEM_NAME],
                    closestNeighborDistance)
    return starLinks
        
#--------------------------------------    
def LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, onlyHab, numNeighbors):
    """call LinkNeighbor_ numNeighbors times"""
    # insert dlog to get numNeighbors
    for x in range(numNeighbors):
        starLinks = LinkNeighbor_(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, onlyHab)
    return starLinks

#--------------------------------------    
def LinkMaxAndMin_(aStarSystem, starSystems, starLinks, linkFile, minDist, 
                    maxDist, maxLinks):
    """given aStarSystem, link to all other systems that are farther than
    minDist and closer than maxDist"""
    # insert dlog to get minDist, maxDist, and maxLinks
    numLinks = 0
    for systemIndex in starSystems:
        destSystem = starSystems[systemIndex]
        distance = TrueDistance(aStarSystem[SystemField.LOCALE], 
                                destSystem[SystemField.LOCALE])
        if (distance >= minDist) and (distance <= maxDist):
            numLinks += 1
            if numLinks <= maxLinks:
                starLinks = AddLink(starLinks, linkFile, 
                        aStarSystem[SystemField.INDEX],
                        aStarSystem[SystemField.SYSTEM_NAME],
                        destSystem[SystemField.INDEX],
                        destSystem[SystemField.SYSTEM_NAME],
                        distance)
            else:
                break
    return starLinks
    
# START OF LINK FUNCTIONS
# these are passable to CalculateLinks as linkFunction

# link to 1 closest neighbor
def LinkNeighbor(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, False, 1)
    return starLinks
    
# link to 1 closest habitable neighbor 
def LinkHabNeighbor(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, True, 1)
    return starLinks
    
# link to 2 closest neighbor
def LinkTwoNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, False, 2)
    return starLinks
    
# link to 2 closest habitable neighbor 
def LinkTwoHabNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, True, 2)
    return starLinks
    
# link to maximum of 6 closest neighbors starting at 5.0 units (parsecs) to a maximum of 10 units
# edit the LinkMaxAndMin_() below to change the values
def LinkMaxAndMin(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkMaxAndMin_(aStarSystem, starSystems, starLinks, linkFile, 5.0, 10.0, 6)
    return starLinks
    
# link to maximum of 100 closest neighbors within 2.5 units (parsecs)
# edit the LinkMaxAndMin_() below to change the values
def LinkDistanceLimit(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkMaxAndMin_(aStarSystem, starSystems, starLinks, linkFile, 0.0, 2.5, 100)
    return starLinks
   
# link to maximum of 100 closest neighbors within 2.5 units, 
# and to closest neighbor regardless of distance
def LinkDistLimitAndNeighbor(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkDistanceLimit(aStarSystem, starSystems, starLinks, linkFile)
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, False, 1)
    return starLinks
    
# link to maximum of 100 closest neighbors within 2.5 units, 
# and to closest 2 neighbor regardless of distance
def LinkDistLimitAndTwoNeighbor(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkDistanceLimit(aStarSystem, starSystems, starLinks, linkFile)
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, False, 2)
    return starLinks
    
# link to maximum of 100 closest neighbors within 2.5 units, 
# and to closest Habitable neighbor regardless of distance
def LinkDistLimitAndHabNeighbor(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkDistanceLimit(aStarSystem, starSystems, starLinks, linkFile)
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, True, 1)
    return starLinks

# link to maximum of 100 closest neighbors within 2.5 units, 
# and to closest 2 Habitable neighbor regardless of distance
def LinkDistLimitAndTwoHab(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkDistanceLimit(aStarSystem, starSystems, starLinks, linkFile)
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, True, 2)
    return starLinks

# link to 2 closest neighbor
# and link to 2 closest habitable neighbor 
def LinkTwoHabAndTwoNeighbor(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars):
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, False, 2)
    starLinks = LinkNNeighbors(aStarSystem, starSystems, starLinks, linkFile, forbiddenStars, True, 2)
    return starLinks
    
# END OF LINK FUNCTIONS

#--------------------------------------    
def CalculateLinks(theFileName, starSystems, linkFunction):
    """this function iterates though the star list, passing
    each star to the selected link function"""
    linkFile = open(theFileName + "Links.txt", 'w')
    starLinks = {}
    index = -1
    numStars = len(starSystems)
    for systemIndex in starSystems:
        forbiddenStars = []
        index += 1
        print "Calculating Star %d of %d" % (index, numStars)        
        starLinks = linkFunction(starSystems[systemIndex], starSystems, 
                    starLinks, linkFile, forbiddenStars)
    linkFile.close()
    return starLinks
    
#--------------------------------------    
def PrintGMLFile(theFileName, starSystems, starLinks):
    """writes a *.gml file suitable for yEd.exe"""
    gmlFile = open(theFileName + ".gml", 'w')
    
    # write header
    gmlFile.write("Creator\t\"yFiles\"\n")
    gmlFile.write("Version\t2.0\n")
    gmlFile.write("graph\n")
    gmlFile.write("[\n")
    gmlFile.write("\thierarchic\t1\n")
    gmlFile.write("\tlabel\t\"\"\n")
    gmlFile.write("\tdirected\t1\n")
    
    # write nodes
    index = -1
    numStars = len(starSystems)
    for systemIndex in starSystems:
        index += 1
        print "Printing Star %d of %d" % (index, numStars)        
        aStarSystem = starSystems[systemIndex]
        nameBuffer = aStarSystem[SystemField.SYSTEM_NAME].replace(
                            ", ", "\n")
        hasAtLeastOneLink = False
        for linkIndex in starLinks:
            aStarLink = starLinks[linkIndex]
            if (aStarLink[LinkField.START_INDEX] == systemIndex or
                aStarLink[LinkField.DEST_INDEX] == systemIndex):
                hasAtLeastOneLink = True
                break
        if hasAtLeastOneLink:
            gmlFile.write("\tnode\n")
            gmlFile.write("\t[\n")
            
            gmlFile.write("\t\tid\t%d\n" % (systemIndex))
            gmlFile.write("\t\tlabel\t\"%s\"\n" % (nameBuffer))
            gmlFile.write("\t\tgraphics\n")
            gmlFile.write("\t\t[\n")

            gmlFile.write("\t\t\tx\t0.0\n")
            gmlFile.write("\t\t\ty\t0.0\n")
            gmlFile.write("\t\t\tw\t35.0\n")
            gmlFile.write("\t\t\th\t35.0\n")

            if aStarSystem[SystemField.AT_LEAST_ONE_HABITABLE]:
                gmlFile.write("\t\t\ttype\t\"ellipse\"\n")
                gmlFile.write("\t\t\twidth\t3.0\n")
                gmlFile.write("\t\t\toutline\t\"#000000\"\n")
            else:
                gmlFile.write("\t\t\ttype\t\"rectangle\"\n")
                gmlFile.write("\t\t\twidth\t1.0\n")
                gmlFile.write("\t\t\toutline\t\"#000000\"\n")

            gmlFile.write("\t\t\tfill\t\"%s\"\n" % (StarNodeColors.color[
                            aStarSystem[SystemField.BRIGHTEST_SPECTRALCLASS]]
                            ))
                            
            gmlFile.write("\t\t]\n")
                
            gmlFile.write("\t\tLabelGraphics\n")
            gmlFile.write("\t\t[\n")

            gmlFile.write("\t\t\ttext\t\"%s\"\n" % (nameBuffer))

            gmlFile.write("\t\t\tfontSize\t12\n")
            gmlFile.write("\t\t\tfontName\t\"Arial Narrow\"\n")
            gmlFile.write("\t\t\tanchor\t\t\"c\"\n")

            gmlFile.write("\t\t]\n")

            gmlFile.write("\t]\n")

    # write edges
    index = -1
    numLinks = len(starLinks)
    for linkIndex in starLinks:
        index += 1
        print "Printing link %d of %d" % (index, numLinks)        
        aStarLink = starLinks[linkIndex]
        startSystem = starSystems[aStarLink[LinkField.START_INDEX]]
        destSystem = starSystems[aStarLink[LinkField.DEST_INDEX]]

        gmlFile.write("\tedge\n")
        gmlFile.write("\t[\n")

        gmlFile.write("\t\tsource\t%d\n" % (aStarLink[LinkField.START_INDEX]))
        gmlFile.write("\t\ttarget\t%d\n" % (aStarLink[LinkField.DEST_INDEX]))

        gmlFile.write("\t\tgraphics\n")
        gmlFile.write("\t\t[\n")

        if (startSystem[SystemField.AT_LEAST_ONE_HABITABLE] and
            destSystem[SystemField.AT_LEAST_ONE_HABITABLE]):
            gmlFile.write("\t\t\twidth\t3\n")
        else:
            gmlFile.write("\t\t\twidth\t1\n")
            
        gmlFile.write("\t\t\ttype\t\"line\"\n")
        gmlFile.write("\t\t\tfill\t\"#000000\"\n")

        gmlFile.write("\t\t]\n")

        gmlFile.write("\t\tLabelGraphics\n")
        gmlFile.write("\t\t[\n")

        gmlFile.write("\t\t\ttext\t\"%.1f\"\n" % 
                        (aStarLink[LinkField.DISTANCE]))

        gmlFile.write("\t\t\toutline\t\"#000000\"\n")
        gmlFile.write("\t\t\tfill\t\t\"#FFFFFF\"\n")
        gmlFile.write("\t\t\tfontSize\t9\n")
        gmlFile.write("\t\t\tfontName\t\"Dialog\"\n")
        gmlFile.write("\t\t\tmodel\t\t\"centered\"\n")
        gmlFile.write("\t\t\tposition\t\"center\"\n")

        gmlFile.write("\t\t]\n")

        gmlFile.write("\t]\n")
            
    # all done
    gmlFile.write("]\n")
    
    gmlFile.close()
    

#******************************************************************************
# module/program switch

if __name__ == "__main__":
    
    parser = optparse.OptionParser("usage: %prog [options] myStarFile.csv")
    parser.add_option("-L", "--link", action="store", type="int", dest="linkCode", default=0, help="star linking algorithm (see help file)")
    (options, args) = parser.parse_args()
    
    if len(args) == 0:
        starfileName = "HabHYG50ly.csv"
        
    starfilePrefix = starfileName[:-4]
        
    linkCode = options.linkCode
    
    linkFunction = LinkNeighbor
    if linkCode == 0:
        linkFunction = LinkNeighbor
    elif linkCode == 1:
        linkFunction = LinkHabNeighbor
    elif linkCode == 2:
        linkFunction = LinkTwoNeighbors
    elif linkCode == 3:
        linkFunction = LinkTwoHabNeighbors
    elif linkCode == 4:
        linkFunction = LinkMaxAndMin
    elif linkCode == 5:
        linkFunction = LinkDistanceLimit
    elif linkCode == 6:
        linkFunction = LinkDistLimitAndNeighbor
    elif linkCode == 7:
        linkFunction = LinkDistLimitAndTwoNeighbor
    elif linkCode == 8:
        linkFunction = LinkDistLimitAndHabNeighbor
    elif linkCode == 9:
        linkFunction = LinkDistLimitAndTwoHab
    elif linkCode == 10:
        linkFunction = LinkTwoHabAndTwoNeighbor
    else:
        linkFunction = LinkNeighbor
        
    print "** Loading Stars"
    starSystems, boundingBox = LoadStars(starfilePrefix)
    print "** Calculating Links"
    starLinks = CalculateLinks(starfilePrefix, starSystems, linkFunction)
    print "** Printing GML File"
    PrintGMLFile(starfilePrefix, starSystems, starLinks)
    print "** all done"
    
