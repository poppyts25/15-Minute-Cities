# import necessary libraries and packages
import warnings
warnings.filterwarnings("ignore")
from flask import Flask, render_template, request, url_for
import osmnx as ox
import networkx as nx
import folium
from folium.plugins import MarkerCluster  
from shapely.geometry import MultiPolygon, Point
import geopandas as gpd
import sys
import numpy as np
from alphashape import alphashape
import os
from requests.exceptions import RequestException


app = Flask(__name__)

def save_map(m, map_path):
    # Get the generated map variable name (Folium uses a unique name)
    map_name = m.get_name()

    # Inject JavaScript to enable map clicks
    click_js = f"""
        document.addEventListener("DOMContentLoaded", function() {{
            var leafletMap = window["{map_name}"];
            if (!leafletMap) {{
                console.error("Map instance not found! Check if Folium initialized the map.");
                return;
            }}

            var currentMarker = null; // Variable to store the current marker

            leafletMap.on('click', function(e) {{
                var coords = e.latlng;
                var lat = coords.lat;
                var lng = coords.lng;

                if (currentMarker) {{
                    leafletMap.removeLayer(currentMarker);
                }}

                currentMarker = L.marker([lat, lng]).addTo(leafletMap);
                currentMarker.bindPopup("Coordinates: " + lat.toFixed(5) + ", " + lng.toFixed(5)).openPopup();

                console.log("Coordinates: [" + lat + ", " + lng + "]");
                // Send coordinates to the parent window (if the map is in an iframe)
                window.parent.postMessage({{ lat: lat, lng: lng }}, "*");
            }});
            console.log("Map click listener added successfully.");
        }});
    """

    m.get_root().html.add_child(folium.Element(f"""
        <script>
        {click_js}
        </script>
    """))

    # Define bounds: [southwest, northeast]
    bounds = [[-60, -180], [60, 180]]  # Approximate bounds around the globe

    m.save(map_path)  # Save in 'static' folder

def add_amenities(m, centre, polygon, amenity, distance):
    if amenity == "park":
        tags = {"leisure": "park"}
    elif amenity == "supermarket":
        tags = {"shop": "supermarket"}
    elif amenity == "bakery":
        tags = {"shop": "bakery"}
    elif amenity == "convenience":
        tags = {"shop": "convenience"}
    else:
        tags = {"amenity": amenity}
    
    epsg_code = get_epsg_code(centre[0], centre[1])
    pois = ox.geometries_from_point(centre, tags=tags, dist=distance)
    isochrone_gdf = gpd.GeoDataFrame({'id': [1], 'geometry': [polygon]}, crs=epsg_code)
    pois = gpd.GeoDataFrame(pois, geometry='geometry', crs=epsg_code)
    pois_filtered = gpd.sjoin(pois, isochrone_gdf, predicate="within", how="inner")
    valid_amenities = pois_filtered[pois_filtered.geometry.notnull()]
    valid_amenities['geometry'] = valid_amenities.geometry.apply(
        lambda geom: geom.centroid if geom.geom_type != 'Point' else geom
    )

    # Add a MarkerCluster layer to the map
    marker_cluster = MarkerCluster().add_to(m)

    # Add each marker to the cluster
    for _, row in valid_amenities.iterrows():
        folium.Marker(
            location=[row.geometry.y, row.geometry.x],
            popup=row.get('name', 'Unnamed Amenity'),
            icon=folium.Icon(color="blue")
        ).add_to(marker_cluster)  # Add to the cluster

    return m

def count_markers(m):
    markers = 0
    map_html = m.get_root().render()

    # Count occurrences of marker instances
    markers = map_html.count('L.marker(')
    
    return markers

def get_epsg_code(lat, lon):
    utm_zone = int((lon + 180) / 6) + 1  # Calculate UTM Zone
    epsg_code = 32600 + utm_zone if lat >= 0 else 32700 + utm_zone  # 326XX for Northern Hemisphere, 327XX for Southern
    return epsg_code

def get_area(isochrone_polygon, centre):
    if isinstance(isochrone_polygon, MultiPolygon):
        combined_polygon = isochrone_polygon
    else:
        combined_polygon = MultiPolygon([isochrone_polygon])
    
    # **Determine the correct UTM zone dynamically**
    epsg_code = get_epsg_code(centre[0], centre[1])
    print("EPSG code:", epsg_code)
    # Convert to GeoDataFrame
    isochrone_gdf = gpd.GeoDataFrame({"geometry": [combined_polygon]}, crs="EPSG:4326")

    # **Reproject to UTM for accurate area calculation**
    isochrone_gdf = isochrone_gdf.to_crs(epsg=epsg_code)
    # **Calculate the area in square metres**
    area = isochrone_gdf.geometry.area.iloc[0]/1000000
    print("Area in square metres:", area)
    area = round(area, 2)
    return area

def add_error_alert(m, error_message):
    # Add an alert to the map
    m.get_root().html.add_child(folium.Element(f"""
        <script>
        alert("{error_message}");
        </script>
    """))

@app.route('/')
def index():
    # Generate the map
    map_path = "static/map.html"
    
    # Ensure static directory exists
    if not os.path.exists("static"):
        os.makedirs("static")

    # Create a Folium map
    m = folium.Map(location=[50.7260, -3.5275], zoom_start=12, max_bounds=True ,min_zoom=2 )  # Example: exeter, UK

    # save the map with the click event
    save_map(m, map_path)

    form_data = {
        "travel_method": None,
        "speed": None,
        "amenity": None,
        "latitude": None,
        "longitude": None
    }

    stats = {
        "markers": None,
        "area": None
        }

    return render_template('index.html', map_url=map_path, form_data=form_data, stats=stats)

# Main route for the index page
@app.route("/", methods=["GET", "POST"])
@app.route("/home", methods=["GET", "POST"])
def main_page():
    map_path = "static/map.html"

    if request.method == "POST":
        # Retrieve user input from the form
        travel_method = request.form.get("travel_method")
        speed = request.form.get("speed")
        amenity = request.form.get("amenity")

        # If latitude or longitude is not provided set to exeter city centre
        try:
            centre_lat = float(request.form.get("latitude"))
            centre_lon = float(request.form.get("longitude"))
        except (ValueError, TypeError):
            centre_lat = 50.7260
            centre_lon = -3.5275  # Example: exeter, UK
        
        centre_point = (centre_lat, centre_lon)
        
        # Create a Folium map
        m = folium.Map(location=centre_point, zoom_start=13, max_bounds=True, min_zoom=2)  


        # Define the speed based on user input
        travel_time = 15  # in minutes
        if travel_method == "Walk":
            travel_speed = 3 if speed == "S" else 5 if speed == "A" else 7
            distance = 3000
            travel = "walk"
        else:
            travel_speed = 14 if speed == "S" else 16 if speed == "A" else 18
            distance =  5000 if speed == "S" else 7000 if speed == "A" else 9000
            travel = "bike"
            
        speed_mps = travel_speed * 1000 / 3600  # Convert to metres per second
        travel_time_seconds = 15 * 60  # 15 minutes
        travel_radius = travel_time_seconds * speed_mps

        ## Load network graph
        ox.settings.use_cache = True  # Enable caching to avoid repeated downloads
        try:
            G = ox.graph_from_point((centre_lat, centre_lon), network_type=travel, dist=distance, simplify=True, dist_type="network")
        except Exception as e:
            print(f"Error loading graph: {e}")           

            error_message = "Error loading the network graph. There are no roads in this area or the area is too small. Please try a different location."
            add_error_alert(m, error_message)
            save_map(m, map_path)
            G = None    

            form_data = {
                "travel_method": request.form.get("travel_method", ""),
                "speed": request.form.get("speed", ""),
                "amenity": request.form.get("amenity", ""),
                "latitude": None,
                "longitude": None
            }

            stats = {
                "markers": None,
                "area": None
                }

            return render_template('index.html', map_url=map_path, form_data=form_data, stats=stats)
            

        if G is not None:
            centre_node = ox.distance.nearest_nodes(G, centre_lon, centre_lat)

            isochrone_subgraph = nx.ego_graph(G, centre_node, radius=travel_radius, distance="length")
            nodes, edges = ox.graph_to_gdfs(isochrone_subgraph, nodes=True, edges=True)
            # Extract node coordinates from the isochrone subgraph
            node_points = np.array(list(zip(nodes.geometry.x, nodes.geometry.y)))

            # Adjust alpha parameter (lower values give more detail)
            alpha = 0.3  # Try tweaking between 0.1 and 1.0 for best results
            isochrone_polygon = alphashape(node_points, alpha)

          
            # Add the isochrone polygon to the map
            style_func = lambda x: {"fillColor": "red", "color": "red", "weight": 0, "fillOpacity": 0.4}
            if isinstance(isochrone_polygon, MultiPolygon):
                for poly in isochrone_polygon.geoms:
                    folium.GeoJson(poly, style_function=style_func).add_to(m)
            else:
                folium.GeoJson(isochrone_polygon, style_function=style_func).add_to(m)


            markers=0
            if amenity != "none":
                m = add_amenities(m, centre_point, isochrone_polygon, amenity, distance)
                markers = count_markers(m)

            # Add marker for centre point, selected by user or default centre
            folium.Marker(
                location=[centre_lat, centre_lon],
                popup="Centre Point",
                icon=folium.Icon(color="lightgray")
            ).add_to(m)


            # Serve the map from the static directory
            # save the map with the click event
            save_map(m, map_path)

    

    form_data = {
        "travel_method": request.form.get("travel_method", ""),
        "speed": request.form.get("speed", ""),
        "amenity": request.form.get("amenity", ""),
        "latitude": request.form.get("latitude", ""),
        "longitude": request.form.get("longitude", "")
    }

    area = get_area(isochrone_polygon, centre_point)

    stats = {
        "markers": markers, 
        "area":area
        }

    # Render the index.html page with the map URL
    return render_template("index.html", map_url=map_path, form_data=form_data, stats=stats)


# Main route for the compare areas page
@app.route("/compare")
def compare_index():
    # Generate the map
    map_path_1 = "static/area1.html"
    map_path_2 =  "static/area2.html"

    # Ensure static directory exists
    if not os.path.exists("static"):
        os.makedirs("static")

    # Create a Folium map
    m1 = folium.Map(location=[50.7260, -3.5275], zoom_start=12,max_bounds=True, min_zoom=2)  # Example: exeter, UK
    m2 = folium.Map(location=[50.7260, -3.5275], zoom_start=12, max_bounds=True, min_zoom=2)  # Example: exeter, UK


    save_map(m1, map_path_1)
    save_map(m2, map_path_2)


    form_data = {
        "travel_method": None,
        "speed": None,
        "amenity": None,
        "latitude": None,
        "longitude": None
    }

    area = 0

    stats = {
        "markers": None,
        "area": None
        }
    
    return render_template('compare_area.html', map_url_1=map_path_1, map_url_2=map_path_2, form_data=form_data, stats=stats)

# compare areas page
@app.route("/compare", methods=["GET", "POST"])
def compare_areas():
    map_path_1 = "static/area1.html"  # Serve the map from the static directory
    map_path_2 = "static/area2.html"
    if request.method == "POST":
        # Retrieve user input from the form
        travel_method = request.form.get("travel_method")
        speed = request.form.get("speed")
        amenity = request.form.get("amenity")
        
        try:
            centre_lat_1 = float(request.form.get("latitude1"))
            centre_lon_1 = float(request.form.get("longitude1"))
        except (ValueError, TypeError):
            centre_lat_1 = 50.7260
            centre_lon_1 = -3.5275

        centre_point_1 = (centre_lat_1, centre_lon_1)

        try:
            centre_lat_2 = float(request.form.get("latitude2"))
            centre_lon_2 = float(request.form.get("longitude2"))
        except (ValueError, TypeError):
            centre_lat_2 = 50.7260
            centre_lon_2 = -3.5275

        centre_point_2 = (centre_lat_2, centre_lon_2)

        # Define the speed based on user input
        travel_time = 15  # in minutes
        if travel_method == "Walk":
            travel_speed = 3 if speed == "S" else 5 if speed == "A" else 7
            distance = 3000
            travel = "walk"
        else:
            travel_speed = 14 if speed == "S" else 16 if speed == "A" else 18
            distance =  5000 if speed == "S" else 7000 if speed == "A" else 9000
            travel = "bike"
        speed_mps = travel_speed * 1000 / 3600  # Convert to metres per second
        travel_time_seconds = 15 * 60  # 15 minutes
        travel_radius = travel_time_seconds * speed_mps


        ## Load network graph
        ox.settings.use_cache = True  # Enable caching to avoid repeated downloads
        
        try:
            G1 = ox.graph_from_point((centre_lat_1, centre_lon_1), network_type=travel, dist=distance, simplify=True, dist_type="network")
        except Exception as e:
            print(f"Error loading graph for area 1: {e}")           

            error_message = "Error loading the network graph for area 1. There are no roads in this area or the area is too small. Please try a different location."
            m1 = folium.Map(location=centre_point_1, zoom_start=13, max_bounds=True, min_zoom=2)
            add_error_alert(m1, error_message)
            save_map(m1, map_path_1)
            G1 = None

            form_data = {
                "travel_method": request.form.get("travel_method", ""),
                "speed": request.form.get("speed", ""),
                "amenity": request.form.get("amenity", ""),
                "latitude1": request.form.get("latitude1", ""),
                "longitude1": request.form.get("longitude1", ""),
                "latitude2": request.form.get("latitude2", ""),
                "longitude2": request.form.get("longitude2", "")
            }


            stats = {
                "markers_1": None,
                "markers_2": None,
                "area_1": None,
                "area_2": None
            }

            # Render the compare.html page with the map URL
            return render_template('compare_area.html', map_url_1=map_path_1, map_url_2=map_path_2, form_data=form_data, stats=stats)
            
        try:
            G2 = ox.graph_from_point((centre_lat_2, centre_lon_2), network_type=travel, dist=distance, simplify=True, dist_type="network")
        except Exception as e:
            print(f"Error loading graph for area 2: {e}")           

            error_message = "Error loading the network graph for area 2. There are no roads in this area or the area is too small. Please try a different location."
            m2 = folium.Map(location=centre_point_2, zoom_start=13, max_bounds=True, min_zoom=2)
            add_error_alert(m2, error_message)
            save_map(m2, map_path_2)
            G2 = None


            form_data = {
                "travel_method": request.form.get("travel_method", ""),
                "speed": request.form.get("speed", ""),
                "amenity": request.form.get("amenity", ""),
                "latitude1": request.form.get("latitude1", ""),
                "longitude1": request.form.get("longitude1", ""),
                "latitude2": request.form.get("latitude2", ""),
                "longitude2": request.form.get("longitude2", "")
            }

            stats = {
                "markers_1": None,
                "markers_2": None,
                "area_1": None,
                "area_2": None
                }

            # Render the compare.html page with the map URL
            return render_template('compare_area.html', map_url_1=map_path_1, map_url_2=map_path_2, form_data=form_data, stats=stats)

        if (G1 is not None) and (G2 is not None):
            centre_node_1 = ox.distance.nearest_nodes(G1, centre_lon_1, centre_lat_1)
            centre_node_2 = ox.distance.nearest_nodes(G2, centre_lon_2, centre_lat_2)

            isochrone_subgraph_1 = nx.ego_graph(G1, centre_node_1, radius=travel_radius, distance="length")
            isochrone_subgraph_2 = nx.ego_graph(G2, centre_node_2, radius=travel_radius, distance="length")

            nodes_1, edges_1 = ox.graph_to_gdfs(isochrone_subgraph_1, nodes=True, edges=True)
            nodes_2, edges_2 = ox.graph_to_gdfs(isochrone_subgraph_2, nodes=True, edges=True)
            # Extract node coordinates from the isochrone subgraph
            node_points_1 = np.array(list(zip(nodes_1.geometry.x, nodes_1.geometry.y)))
            node_points_2 = np.array(list(zip(nodes_2.geometry.x, nodes_2.geometry.y)))

            # Adjust alpha parameter (lower values give more detail)
            alpha = 0.3 
            isochrone_polygon_1 = alphashape(node_points_1, alpha)
            isochrone_polygon_2 = alphashape(node_points_2, alpha)

            # Create a Folium map
            m1 = folium.Map(location=centre_point_1, zoom_start=13,max_bounds=True, min_zoom=2)
            m2 = folium.Map(location=centre_point_2, zoom_start=13, max_bounds=True, min_zoom=2 )

            # Add the isochrone polygon to the map
            style_func = lambda x: {"fillColor": "red", "color": "red", "weight": 0, "fillOpacity": 0.4}
            if isinstance(isochrone_polygon_1, MultiPolygon):
                for poly in isochrone_polygon_1.geoms:
                    folium.GeoJson(poly, style_function=style_func).add_to(m1)
            else:
                folium.GeoJson(isochrone_polygon_1, style_function=style_func).add_to(m1)

            if isinstance(isochrone_polygon_2, MultiPolygon):
                for poly in isochrone_polygon_2.geoms:
                    folium.GeoJson(poly, style_function=style_func).add_to(m2)
            else:
                folium.GeoJson(isochrone_polygon_2, style_function=style_func).add_to(m2)


            markers_1=0
            markers_2=0
            # (Optional) Add amenity markers if applicable
            if amenity != "none":
                m1 = add_amenities(m1, centre_point_1, isochrone_polygon_1, amenity, distance)
                markers_1 = count_markers(m1)

                m2 = add_amenities(m2, centre_point_2, isochrone_polygon_2, amenity, distance)
                markers_2 = count_markers(m2)


            # Add marker for centre point, selected by user or default centre
            folium.Marker(
                location=[centre_lat_1, centre_lon_1],
                popup="Centre Point",
                icon=folium.Icon(color="lightgray")
            ).add_to(m1)


            folium.Marker(
                location=[centre_lat_2, centre_lon_2],
                popup="Centre Point",
                icon=folium.Icon(color="lightgray")
            ).add_to(m2)

            
            save_map(m1, map_path_1)
            save_map(m2, map_path_2)


    form_data = {
        "travel_method": request.form.get("travel_method", ""),
        "speed": request.form.get("speed", ""),
        "amenity": request.form.get("amenity", ""),
        "latitude1": request.form.get("latitude1", ""),
        "longitude1": request.form.get("longitude1", ""),
        "latitude2": request.form.get("latitude2", ""),
        "longitude2": request.form.get("longitude2", "")
    }

    area_1 = get_area(isochrone_polygon_1, centre_point_1)
    area_2 = get_area(isochrone_polygon_2, centre_point_2)

    stats = {
        "markers_1": markers_1,
        "markers_2": markers_2,
        "area_1": area_1,
        "area_2": area_2
        }

    # Render the compare.html page with the map URL
    return render_template('compare_area.html', map_url_1=map_path_1, map_url_2=map_path_2, form_data=form_data, stats=stats)
    
# about page
@app.route("/about")
def about():
    return render_template('about.html')

if __name__ == "__main__":
    warnings.simplefilter("ignore")
    app.run(debug=True, threaded=True)
