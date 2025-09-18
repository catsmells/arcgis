import arcpy
from queue import PriorityQueue
import time
import math
class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the .pyt file)."""
        self.label = "Pole of Inaccessibility Toolbox"
        self.alias = "poi"
        self.tools = [PoleOfInaccessibilityTool]
class PoleOfInaccessibilityTool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Calculate Pole of Inaccessibility"
        self.description = "Calculates the pole of inaccessibility for a given polygon area."
        self.canRunInBackground = False
    def getParameterInfo(self):
        """Define parameter definitions"""
        input_polygon = arcpy.Parameter(
            displayName="Input Polygon",
            name="input_polygon",
            datatype="FeatureSet",
            parameterType="Required",
            direction="Input"
        )
        # To allow drawing polygons, note: In ArcGIS Pro, FeatureSet will allow interactive drawing.
        # For polygon type, you may need to set a schema from a template polygon feature class.
        # Example: Create an empty polygon FC and set input_polygon.schema = arcpy.Describe("path_to_template_fc").shapeFieldName or similar.
        # For simplicity, assuming user draws polygon.
        output_point = arcpy.Parameter(
            displayName="Output Point Feature Class",
            name="output_point",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        tolerance = arcpy.Parameter(
            displayName="Tolerance",
            name="tolerance",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )
        tolerance.value = 0.001  # Default tolerance

        return [input_polygon, output_point, tolerance]

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal validation is performed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool parameter."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        input_features = parameters[0].value
        tolerance = parameters[2].value if parameters[2].value else 0.001

        # Copy input features to in-memory feature class
        temp_fc = arcpy.management.CopyFeatures(input_features, "in_memory/temp_polygon")

        # Assume single polygon for simplicity; get the geometry
        with arcpy.da.SearchCursor(temp_fc, ["SHAPE@"]) as cursor:
            for row in cursor:
                polygon_geom = row[0]
                break  # Take first feature

        # Compute bounding box
        extent = polygon_geom.extent
        min_x = extent.XMin
        min_y = extent.YMin
        max_x = extent.XMax
        max_y = extent.YMax
        width = max_x - min_x
        height = max_y - min_y
        cell_size = min(width, height)
        h = cell_size / 2.0
        max_dim = max(width, height)

        if cell_size == 0 or polygon_geom.area < max_dim * tolerance:
            # Degenerate case
            pole_pt = polygon_geom.centroid
        else:
            # Define Cell class
            class Cell:
                def __init__(self, x, y, h, polygon_geom):
                    self.x = x
                    self.y = y
                    self.h = h
                    pt = arcpy.Point(x, y)
                    pt_geom = arcpy.PointGeometry(pt, polygon_geom.spatialReference)
                    self.d = pt_geom.distanceTo(polygon_geom)  # Distance to boundary
                    self.max = self.d + self.h * math.sqrt(2)

            # Priority queue for cells (max heap using negative max)
            cell_queue = PriorityQueue()

            # Cover polygon with initial cells
            x = min_x
            while x < max_x:
                y = min_y
                while y < max_y:
                    c = Cell(x + h, y + h, h, polygon_geom)
                    cell_queue.put((-c.max, time.time(), c))
                    y += cell_size
                x += cell_size

            # Initial best cell at centroid
            def get_centroid_cell(polygon_geom):
                cent_pt = polygon_geom.centroid
                return Cell(cent_pt.X, cent_pt.Y, 0, polygon_geom)

            best_cell = get_centroid_cell(polygon_geom)

            # Check bbox center
            bbox_cell = Cell(min_x + width / 2, min_y + height / 2, 0, polygon_geom)
            if bbox_cell.d > best_cell.d:
                best_cell = bbox_cell

            # Iterate
            while not cell_queue.empty():
                _, _, cell = cell_queue.get()

                if cell.d > best_cell.d:
                    best_cell = cell

                if cell.max - best_cell.d <= tolerance:
                    continue

                h = cell.h / 2
                for dx, dy in [(-h, -h), (h, -h), (-h, h), (h, h)]:
                    c = Cell(cell.x + dx, cell.y + dy, h, polygon_geom)
                    cell_queue.put((-c.max, time.time(), c))

            pole_pt = arcpy.Point(best_cell.x, best_cell.y)
        # Create output point feature class
        sr = polygon_geom.spatialReference
        output_fc = arcpy.management.CreateFeatureclass("in_memory", "output_point", "POINT", spatial_reference=sr)
        with arcpy.da.InsertCursor(output_fc, ["SHAPE@"]) as cursor:
            cursor.insertRow([pole_pt])
        # Set the derived output parameter
        parameters[1].value = output_fc
        arcpy.management.Delete("in_memory")
        return
