# Test suite to test functions in app.py

import pytest
import folium
import os
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch
from shapely.geometry import Point, Polygon, MultiPolygon
import geopandas as gpd
from folium.plugins import MarkerCluster

# Import the functions to be tested from app.py
from app import save_map, add_amenities, count_markers, get_epsg_code, get_area

# Test cases for the save_map function in app.py
def test_save_map_creates_file_and_contains_js():
    with TemporaryDirectory() as tmpdir:
        map_path = os.path.join(tmpdir, "test_map.html")
        m = folium.Map(location=[51.5, -0.1], zoom_start=10)

        save_map(m, map_path)

        # Check file is created
        assert os.path.exists(map_path), "Map HTML file was not created."

        # Check file is not empty
        assert os.path.getsize(map_path) > 0, "Map HTML file is empty."

        # Check that custom JS was injected
        with open(map_path, "r", encoding="utf-8") as f:
            content = f.read()
            assert "Map click listener added successfully" in content, "Custom JS not found in map HTML."

def test_save_map_does_not_raise_exceptions():
    with TemporaryDirectory() as tmpdir:
        map_path = os.path.join(tmpdir, "test_map.html")
        m = folium.Map(location=[0, 0], zoom_start=2)

        try:
            save_map(m, map_path)
        except Exception as e:
            pytest.fail(f"save_map raised an unexpected exception: {e}")

# Test case for the add_amenities function
@patch("app.get_epsg_code", return_value="EPSG:4326")
@patch("app.ox.geometries_from_point")
def test_add_amenities_adds_markers(mock_geometries_from_point, mock_get_epsg_code):
    # Setup a simple folium map
    m = folium.Map(location=[51.5, -0.1], zoom_start=12)
    
    # Dummy centre and polygon
    centre = (51.5, -0.1)
    polygon = Polygon([(-0.2, 51.4), (-0.2, 51.6), (0, 51.6), (0, 51.4)])
    distance = 1000

    # Create dummy amenities GeoDataFrame with 2 points
    dummy_data = {
        'geometry': [Point(-0.1, 51.5), Point(-0.12, 51.52)],
        'name': ['Test Amenity 1', 'Test Amenity 2']
    }
    gdf = gpd.GeoDataFrame(dummy_data, crs="EPSG:4326")

    mock_geometries_from_point.return_value = gdf

    # Call the function
    result_map = add_amenities(m, centre, polygon, "park", distance)

    # Assertions
    assert isinstance(result_map, folium.Map), "Returned object is not a folium.Map."

    # Check that MarkerCluster layer is added
    marker_clusters = [
        child for child in result_map._children.values()
        if isinstance(child, MarkerCluster)
    ]
    assert len(marker_clusters) == 1, "MarkerCluster was not added to the map."

    # Check that at least 2 markers are present inside the MarkerCluster
    marker_cluster = marker_clusters[0]
    markers = [
        child for child in marker_cluster._children.values()
        if isinstance(child, folium.map.Marker)
    ]
    assert len(markers) == 2, f"Expected 2 markers, found {len(markers)}."

    # Check marker locations
    marker_locations = [
        marker.location for marker in markers
    ]
    assert any(abs(lat - 51.5) < 0.001 for lat, lon in marker_locations), "Marker near expected location missing."


def test_count_markers_counts_correctly():

    # Create a folium map
    m = folium.Map(location=[51.5, -0.1], zoom_start=13)

    # Initially, there should be 0 markers
    count = count_markers(m)
    assert count == 0, f"Expected 0 markers, found {count}"

    # Add 2 markers
    folium.Marker(location=[51.5, -0.1]).add_to(m)
    folium.Marker(location=[51.51, -0.11]).add_to(m)

    # Now, there should be 2 markers
    count = count_markers(m)
    assert count == 2, f"Expected 2 markers, found {count}"

    # Add another marker
    folium.Marker(location=[51.49, -0.09]).add_to(m)

    # Now, there should be 3 markers
    count = count_markers(m)
    assert count == 3, f"Expected 3 markers, found {count}"


@pytest.mark.parametrize("lat, lon, expected_epsg", [
    (0, 0, 32631),             # Equator, Greenwich, should be UTM zone 31N
    (45, -123, 32610),         # Portland, USA, should be UTM zone 10N
    (-33.9, 151.2, 32756),     # Sydney, Australia, should be UTM zone 56S
    (51.5, -0.1, 32630),       # London, UK, should be UTM zone 30N
    (-23.5, -46.6, 32723),     # São Paulo, Brazil, should be UTM zone 23S
])
def test_get_epsg_code_correctness(lat, lon, expected_epsg):
    epsg = get_epsg_code(lat, lon)
    assert epsg == expected_epsg, f"For lat={lat}, lon={lon}, expected EPSG {expected_epsg}, got {epsg}"

def test_get_area_with_polygon():
    polygon = Polygon([
        (-0.1, 51.5),
        (-0.1, 51.51),
        (-0.09, 51.51),
        (-0.09, 51.5),
        (-0.1, 51.5)
    ])
    centre = (51.505, -0.095)

    area = get_area(polygon, centre)
    # Expected ~0.77 km²
    assert 0.7 < area < 0.85, f"Expected area around 0.77 km², got {area} km²"

def test_get_area_with_multipolygon():
    poly1 = Polygon([
        (-0.1, 51.5),
        (-0.1, 51.51),
        (-0.09, 51.51),
        (-0.09, 51.5),
        (-0.1, 51.5)
    ])
    poly2 = Polygon([
        (-0.09, 51.5),
        (-0.09, 51.51),
        (-0.08, 51.51),
        (-0.08, 51.5),
        (-0.09, 51.5)
    ])
    multipolygon = MultiPolygon([poly1, poly2])
    centre = (51.505, -0.095)

    area = get_area(multipolygon, centre)
    # Expected ~1.54 km²
    assert 1.4 < area < 1.7, f"Expected area around 1.54 km², got {area} km²"


