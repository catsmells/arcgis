from __future__ import annotations

import arcpy
from bisect import insort
from dataclasses import dataclass
from collections.abc import Callable, Iterator

# Iterator needed so multipart polygons can get a POI per part
POIMethod = Callable[[arcpy.Polygon, float], Iterator[arcpy.PointGeometry]]
sqrt_2 = 2**0.5

@dataclass
class Cell:
    x: float
    y: float
    h: float
    polygon: arcpy.Polygon
    boundary: arcpy.Polyline
    
    def __post_init__(self) -> None:
        self.centroid = arcpy.PointGeometry(arcpy.Point(self.x, self.y), self.polygon.spatialReference)
        self.distance = self.centroid.distanceTo(self.boundary)
        self.possible_dist = self.distance + self.h * sqrt_2
        # Negative distance for cells that are outside polygon
        if not self.polygon.contains(self.centroid):
            self.distance = -self.distance
            self.possible_dist = -self.possible_dist
    
    def next_cells(self) -> list[Cell]:
        h = self.h / 2
        x_plus = self.x+h
        y_plus = self.y+h
        y_minus = self.y-h
        x_minus = self.x-h
        
        # Subdivide cell and return a list sorted by longest first
        next_cells: list[Cell] = [
            Cell(x_plus, y_plus, h, self.polygon, self.boundary),
            Cell(x_plus, y_minus, h, self.polygon, self.boundary),
            Cell(x_minus, y_plus, h, self.polygon, self.boundary),
            Cell(x_minus , y_minus, h, self.polygon, self.boundary),
        ]
        return sorted([c for c in next_cells if c])

    def __lt__(self, other: Cell) -> bool:
        return self.possible_dist < other.possible_dist
    
    def __gt__(self, other: Cell) -> bool:
        return self.possible_dist > other.possible_dist
    
    def _bool__(self) -> bool:
        return self.possible_dist > 0
    
def adaptive_grid(polygon: arcpy.Polygon, precision: float=0.1) -> Iterator[arcpy.PointGeometry]:
    # Handle multipart
    if polygon.isMultipart:
        polygons = [arcpy.Polygon(part, polygon.spatialReference) for part in polygon]
        yield from (poi for polygon in polygons for poi in adaptive_grid(polygon, precision))
        return

    # Get initial Cell
    _extent = polygon.extent
    _centroid = polygon.centroid
    _cell_size = min(_extent.width, _extent.height) / 2
    
    base_cell = Cell(_centroid.X, _centroid.Y, _cell_size, polygon, polygon.boundary())
    extent_cell = Cell(_extent.polygon.centroid.X, _extent.polygon.centroid.Y, _cell_size, polygon, polygon.boundary())
    
    # Initialize Queue
    cell_queue: list[Cell] = []
    # Polygon and Extent cells
    insort(cell_queue, base_cell)
    insort(cell_queue, extent_cell)
    
    # Initial cells
    _cell_size /= 2
    for x in range(int(_extent.XMin), int(_extent.XMax), int(_cell_size)):
        for y in range(int(_extent.YMin), int(_extent.YMax), int(_cell_size)):
            c = Cell(x + _cell_size, y + _cell_size, _cell_size, polygon, polygon.boundary())
            if c: # Skip cells that are outside polygon
                insort(cell_queue, c)
    
    # Pop the best cell candidate off the queue
    best_cell = cell_queue.pop()
    
    # Traverse the cell subdivisions until a minimum is found for precision level
    iterations = 0
    while cell_queue:
        arcpy.SetProgressorLabel(f'Best Distance: {best_cell.distance:0.0f} {polygon.spatialReference.linearUnitName} (iteration {iterations})...')
        iterations += 1
        # Pop the next cell off the queue
        next_cell = cell_queue.pop()
        if next_cell.possible_dist < best_cell.distance:
            continue # Skip cells with no possible better distance
        
        # Update best cell if next cell is further away
        if next_cell.distance > best_cell.distance:
            best_cell = next_cell
        
        # Add next cell children to the queue
        for cell in next_cell.next_cells():
            if cell.h >= precision and cell.possible_dist > best_cell.distance:
                insort(cell_queue, cell)
        
    arcpy.AddMessage(f'Best Distance: {best_cell.distance:0.0f} {polygon.spatialReference.linearUnitName} (found in {iterations} iterations)')
    yield best_cell.centroid

def b9_hillclimbing(polygon: arcpy.Polygon, precision: float=0.1) -> Iterator[arcpy.PointGeometry]:
    raise NotImplementedError

METHODS: dict[str, POIMethod] = {
    'Adaptive Grid' : adaptive_grid,
    'B9-Hillclimbing': b9_hillclimbing,
}

class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the .pyt file)."""
        self.label = "Pole of Inaccessibility Toolbox"
        self.alias = "poi_toolbox"
        self.tools = [PoleOfInaccessibilityTool]
        
class PoleOfInaccessibilityTool:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Calculate Pole of Inaccessibility"
        self.description = "Calculates the pole of inaccessibility for a given polygon area."
        self.canRunInBackground = False
        
    def getParameterInfo(self):
        """Define parameter definitions"""
        
        polygon_layer = arcpy.Parameter(
            displayName="Polygon Layer",
            name="polygon_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input",
        )
        # Enforce Polygons
        if polygon_layer.filter:
            polygon_layer.filter.list = ['Polygon']
        
        # Allow feature creation using GPFeatureLayer type (Polygon filter applied)
        polygon_layer.controlCLSID = '{60061247-BCA8-473E-A7AF-A2026DDE1C2D}'

        output_poles = arcpy.Parameter(
            displayName="Output Poles",
            name="output_poles",
            datatype="GPString",
            parameterType="Required",
            direction="Output",
        )
        output_poles.value = 'Poles_Of_Inaccessibility'
        
        precision = arcpy.Parameter(
            displayName="Precision",
            name="precision",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
        )
        precision.value = 0.1  # Default tolerance

        method = arcpy.Parameter(
            displayName='Caclulation Method',
            name='method',
            datatype='GPString',
            parameterType='Required',
            direction='Input',
        )
        if method.filter:
            method.filter.list = list(METHODS.keys())
        method.value = 'Adaptive Grid'

        return [polygon_layer, output_poles, precision, method]

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal validation is performed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool parameter."""
        return

    def execute(self, parameters: list[arcpy.Parameter], messages):
        """The source code of the tool.""" 
        params = {p.name: p for p in parameters}
        
        # Get Params
        input_features = params['polygon_layer'].value
        output_poles = params['output_poles'].valueAsText
        precision = params['precision'].value
        get_poi = METHODS[params['method'].value]
        
        # Get POIs
        arcpy.management.CopyFeatures([ # type: ignore
            poi
            for polygon, in arcpy.da.SearchCursor(input_features, ['SHAPE@'])
            if isinstance(polygon, arcpy.Polygon)
            for poi in get_poi(polygon, precision)
        ], output_poles)
        