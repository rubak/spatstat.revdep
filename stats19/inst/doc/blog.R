## ---- eval=FALSE, message=FALSE-----------------------------------------------
#  # release version - currently 0.2.0
#  install.packages("stats19")
#  # dev version
#  # remotes::install_github("ropensci/stats19")

## -----------------------------------------------------------------------------
library(stats19)

## ----dl2017-accidents, message=FALSE------------------------------------------
crashes_2017 = get_stats19(year = 2017, type = "Accidents", ask = FALSE)
nrow(crashes_2017)

## ----crashes_2017-explore-----------------------------------------------------
column_names = names(crashes_2017)
length(column_names)
head(column_names)
class(crashes_2017)
kableExtra::kable(head(crashes_2017[, c(1, 4, 5, 7, 10)]))

## ----file_names---------------------------------------------------------------
stats19::file_names$dftRoadSafetyData_Vehicles_2017.zip

## -----------------------------------------------------------------------------
crashes_2017_raw = get_stats19(year = 2017, type = "Accidents", ask = FALSE, format = FALSE)

## ----crashes_2017-raw---------------------------------------------------------
kableExtra::kable(cbind(head(crashes_2017_raw[1:2, c(7, 10)]), head(crashes_2017[1:2, c(7, 10)])))
class(crashes_2017_raw$Date)
class(crashes_2017$date)

## -----------------------------------------------------------------------------
class(crashes_2017$date)
class(crashes_2017_raw$Date)

## ----format-crashes-sf--------------------------------------------------------
crashes_sf = format_sf(crashes_2017)
# crashes_sf = format_sf(crashes_2017, lonlat = TRUE) # provides the data in lon/lat format

## ----nfatalities, message=FALSE-----------------------------------------------
library(sf)
library(dplyr)
crashes_sf %>% 
  filter(accident_severity == "Fatal") %>% 
  select(n_fatalities = accident_index) %>% 
  aggregate(by = police_boundaries, FUN = length) %>% 
  plot()

## ----ukboundaries-------------------------------------------------------------
west_yorkshire =
  police_boundaries[police_boundaries$pfa16nm == "West Yorkshire", ]

## ----crashes-west_yorkshire---------------------------------------------------
crashes_wy = crashes_sf[west_yorkshire, ]
nrow(crashes_wy) # which is 3.36%

## ----dl2017-vehcas, message=FALSE---------------------------------------------
#crashes_2017 = get_stats19(year = 2017, type = "Accidents", ask = FALSE)
casualties_2017 = get_stats19(year = 2017, type = "Casualties", ask = FALSE)
nrow(casualties_2017)
vehicles_2017 = get_stats19(year = 2017, type = "Vehicles", ask = FALSE)
nrow(vehicles_2017)

## ----table-join, message = FALSE----------------------------------------------
library(tidyr)
library(dplyr)
sel = casualties_2017$accident_index %in% crashes_wy$accident_index
casualties_wy = casualties_2017[sel, ]
cas_types = casualties_wy %>% 
  select(accident_index, casualty_type) %>% 
  group_by(accident_index) %>% 
  summarise(
    Total = n(),
    walking = sum(casualty_type == "Pedestrian"),
    cycling = sum(casualty_type == "Cyclist"),
    passenger = sum(casualty_type == "Car occupant")
    ) 
cj = left_join(crashes_wy, cas_types)

## ----table-join-examples------------------------------------------------------
base::setdiff(names(cj), names(crashes_wy))

## ----sfplot, fig.show='hold', out.width="100%", fig.cap="Spatial distribution of serious and fatal collisions in which people who were walking on the road network ('pedestrians') were hit by a car or other vehicle.", fig.width=9, fig.height=7----
library(ggplot2)
crashes_types = cj %>% 
  filter(accident_severity != "Slight") %>% 
  mutate(type = case_when(
    walking > 0 ~ "Walking",
    cycling > 0 ~ "Cycling",
    passenger > 0 ~ "Passenger",
    TRUE ~ "Other"
  ))
ggplot(crashes_types, aes(size = Total, colour = speed_limit)) +
  geom_sf(show.legend = "point", alpha = 0.3) +
  facet_grid(vars(type), vars(accident_severity)) +
  scale_size(
    breaks = c(1:3, 12),
    labels = c(1:2, "3+", 12)
    ) +
  scale_color_gradientn(colours = c("blue", "yellow", "red")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## ----crashes-map, fig.show='hold', out.width="100%", fig.cap="Spatial distribution of all collisions in which people who were walking on the road network ('pedestrians') were hit by a car or other vehicle in 2017 within West Yorkshire boundary.", fig.width=9, fig.height=7----
library(leaflet)
crashes_pedestrians = crashes_types %>% 
  filter(walking > 0)
# convert to lon lat CRS
crashes_pedestrians_lonlat = st_transform(crashes_pedestrians, crs = 4326)
pal = colorFactor(palette = "Reds", domain = crashes_pedestrians_lonlat$accident_severity, reverse = TRUE)
map = leaflet(data = crashes_pedestrians_lonlat, height = "280px") %>%
  addProviderTiles(provider = providers$OpenStreetMap.BlackAndWhite) %>%
  addCircleMarkers(radius = 0.5, color = ~pal(accident_severity)) %>% 
  addLegend(pal = pal, values = ~accident_severity) %>% 
  leaflet::addMiniMap(toggleDisplay = TRUE)
# map # if you like to see the leaflet version

## ----custom-leaflet-----------------------------------------------------------
library(geojsonsf)
library(htmltools)
geojson = sf_geojson(
  crashes_pedestrians_lonlat[,c("accident_severity")])
template = paste0('
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.4.0/dist/leaflet.css" integrity="sha512-puBpdR0798OZvTTbP4A8Ix/l+A4dHDD0DGqYW6RQ+9jxkRFclaxxQb/SJAWZfWAkuyeQUytO7+7N4QKrDh+drA==" crossorigin=""/>
<script src="https://unpkg.com/leaflet@1.4.0/dist/leaflet.js" integrity="sha512-QVftwZFqvtRNi0ZyCtsznlKSWOStnDORoefr1enyq5mVL4tmKB3S/EnC3rRJcxCPavG10IcrVGSmPh6Qw5lwrg==" crossorigin=""></script>
<div id="mapid" style="width: 100%; height: 400px;">
<script>
	var map = L.map("mapid");
	L.tileLayer("https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoibWFwYm94IiwiYSI6ImNpejY4NXVycTA2emYycXBndHRqcmZ3N3gifQ.rJcFIG214AriISLbB6B5aw", {
		maxZoom: 18,
		attribution: \'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, <a href="https://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, Imagery © <a href="https://www.mapbox.com/">Mapbox</a>\',
		id: "mapbox.streets"
  }).addTo(map);   
  var json = ', geojson, ';', 
  '
  var geojsonMarkerOptions = {
    radius: 6,
    color: "#000",
    weight: 1,
    opacity: 1,
    fillOpacity: 0.5
  };
  var layer = L.geoJSON(json, {
    pointToLayer: function (feature, latlng) {
        return L.circleMarker(latlng, geojsonMarkerOptions);
    },
    style: function(feature) {
      switch (feature.properties.accident_severity) {
          case "Serious": return {color: "#FEB24C"};
          case "Fatal":   return {color: "#BD0026"};
       }
    }
  }).addTo(map);
  map.fitBounds(layer.getBounds());
  var legend = L.control({position: "bottomright"});
	legend.onAdd = function (map) {
		var div = L.DomUtil.create("div", "info legend"), labels = [];
    labels.push("<i style=\'background:#FEB24C\'></i>Serious");
    labels.push("<i style=\'background:#BD0026\'></i>Fatal");
		div.innerHTML = labels.join("<br>");
		return div;
	};
  legend.addTo(map);
  // control that shows state info on hover
	var info = L.control();
	info.onAdd = function (map) {
		this._div = L.DomUtil.create("div", "info");
		this.update();
		return this._div;
	};
	info.update = function (props) {
		this._div.innerHTML = "<h6>Crashes in West Yorkshire (2017)</h6>";
	};
	info.addTo(map);
</script>
<style>
.info { padding: 6px 8px; font: 14px/16px Arial, Helvetica, sans-serif; background: white; background: rgba(255,255,255,0.8); box-shadow: 0 0 15px rgba(0,0,0,0.2); border-radius: 5px; } .info h4 { margin: 0 0 5px; color: #777; }
.legend { text-align: left; line-height: 18px; color: #555; } .legend i { width: 18px; height: 18px; float: left; margin-right: 8px; opacity: 0.7; }</style>
</div>')
path = file.path(tempdir(), "temp.html")
write(template, path)
includeHTML(path)

