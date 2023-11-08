#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DD4hepDetectorHelper.h"
#include <array>
#include <map>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det    = e;
  Material                     vacuum   = description.vacuum();
  int                          det_id   = x_det.id();
  string                       det_name = x_det.nameStr();
  bool                         reflect  = x_det.reflect(false);
  DetElement                   sdet(det_name, det_id);
  Assembly                     assembly(det_name);

  Material air = description.material("Air");
  Volume                             motherVol = description.pickMotherVolume(sdet);
  int                                m_id = 0, c_id = 0, n_sensor = 0;
  map<string, Volume>                modules;
  map<string, Placements>            sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;
  PlacedVolume                       pv;

  double RSU_length = 21.666 * dd4hep::mm;
  //double RSU_length = 30 * dd4hep::mm;
  double RSU_width = 19.564 * dd4hep::mm;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto &params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(
      sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
						    "boundary_material");
  }

  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");


  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
    xml_comp_t x_mod = mi;
    string     m_nam = x_mod.nameStr();
    //xml_comp_t trd   = x_mod.trd();

    double     posY;
    // double     x1              = trd.x1();
    // double     x2              = trd.x2()*0.9;
    double     x1              = RSU_width / 2;
    double     x2              = RSU_width / 2;
    double     z               = RSU_length / 2;
    double     total_thickness = 0.;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci)
      total_thickness += xml_comp_t(ci).thickness();

    double    thickness_so_far = 0.0;
    double    y1               = total_thickness / 2;
    double    y2               = total_thickness / 2;
    Trapezoid m_solid(x1, x2, y1, y2, z);
    Volume    m_volume(m_nam, m_solid, vacuum);
    m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));

    for (ci.reset(), n_sensor = 1, c_id = 0, posY = -y1; ci; ++ci, ++c_id) {
      xml_comp_t c           = ci;
      double     c_thick     = c.thickness();
      auto       comp_x1     = getAttrOrDefault(c, _Unicode(x1), x1);
      auto       comp_x2     = getAttrOrDefault(c, _Unicode(x2), x2);
      auto       comp_height = getAttrOrDefault(c, _Unicode(height), z);

      Material c_mat  = description.material(c.materialStr());
      string   c_name = _toString(c_id, "component%d");

      Trapezoid comp_s1(comp_x1, comp_x2, c_thick / 2e0, c_thick / 2e0, comp_height);
      Solid     comp_shape = comp_s1;
      Volume c_vol(c_name, comp_shape, c_mat);

      c_vol.setVisAttributes(description.visAttributes(c.visStr()));
      pv = m_volume.placeVolume(c_vol, Position(0, posY + c_thick / 2, 0));
      if (c.isSensitive()) {
        module_thicknesses[m_nam] = {thickness_so_far + c_thick / 2.0,
                                     total_thickness - thickness_so_far - c_thick / 2.0};
        // std::cout << " adding sensitive volume" << c_name << "\n";
        sdet.check(n_sensor > 2, "SiTrackerEndcap2::fromCompact: " + c_name + " Max of 2 modules allowed!");
        pv.addPhysVolID("sensor", n_sensor);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        ++n_sensor;
        // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
        Vector3D u(0., 0., -1.);
        Vector3D v(-1., 0., 0.);
        Vector3D n(0., 1., 0.);
        // Vector3D o( 0. , 0. , 0. ) ;

        // compute the inner and outer thicknesses that need to be assigned to the tracking surface
        // depending on wether the support is above or below the sensor
        double inner_thickness = module_thicknesses[m_nam][0];
        double outer_thickness = module_thicknesses[m_nam][1];

        SurfaceType type(SurfaceType::Sensitive);

        // if( isStripDetector )
        //  type.setProperty( SurfaceType::Measurement1D , true ) ;

        VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
        volplane_surfaces[m_nam].push_back(surf);

        //--------------------------------------------
      }
      posY += c_thick;
      thickness_so_far += c_thick;
    }
    modules[m_nam] = m_volume;
  }

  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer(li);
    int        l_id    = x_layer.id();
    int        mod_num = 1;

    xml_comp_t l_env      = x_layer.child(_U(envelope));
    string     layer_name = det_name + std::string("_layer") + std::to_string(l_id);

    std::string layer_vis      = l_env.attr<std::string>(_Unicode(vis));
    double      layer_rmin     = l_env.attr<double>(_Unicode(rmin));
    double      layer_rmax     = l_env.attr<double>(_Unicode(rmax));
    double      layer_x_offset = getAttrOrDefault(l_env, _U(x_offset), 0.0);
    double      layer_length   = l_env.attr<double>(_Unicode(length));
    double      layer_zstart   = l_env.attr<double>(_Unicode(zstart));
    double      layer_center_z = layer_zstart + layer_length / 2.0;
    // printout(INFO,"ROOTGDMLParse","+++ Read geometry from GDML file file:%s",input.c_str());
    // std::cout << "SiTracker Endcap layer " << l_id << " zstart = " << layer_zstart/dd4hep::mm << "mm ( " <<
    // layer_length/dd4hep::mm << " mm thick )\n";

    std::cout << "SiTracker Endcap layer " << l_id << " x_offset = " << layer_x_offset/dd4hep::mm << "mm ( " << std::endl;

    // Assembly    layer_assembly(layer_name);
    // assembly.placeVolume(layer_assembly);

    // Tube   layer_tub(layer_rmin, layer_rmax, layer_length / 2);
    // Tube   layer_tub(0, layer_rmax, layer_length / 2);
    // Volume layer_vol(layer_name, layer_tub, air); // Create the layer envelope volume.
    // layer_vol.setVisAttributes(description.visAttributes(layer_vis));

    Tube   layer_tub(0., layer_rmax, layer_length / 2);
    Tube   layer_tub_cutout(0., layer_rmin, layer_length / 2);
    // x axis in geoviewer is flipped? multiply by -1 so disk envelope is snug around the beampipe
    SubtractionSolid layer_shape(layer_tub, layer_tub_cutout, Transform3D(RotationZYX(0,0,0), Position(layer_x_offset,0,0)));
    Volume layer_vol(layer_name, layer_shape, air); // Create the layer envelope volume.
    layer_vol.setVisAttributes(description.visAttributes(layer_vis));

    PlacedVolume layer_pv;
    if (reflect) {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, -layer_center_z));
      // layer_pv =
      // assembly.placeVolume(layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -layer_center_z)));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_N";
    } else {
      // layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, layer_center_z));
      layer_pv =
	assembly.placeVolume(layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, layer_center_z)));
      layer_pv.addPhysVolID("layer", l_id);
      layer_name += "_P";
    }
    DetElement layer_element(sdet, layer_name, l_id);
    layer_element.setPlacement(layer_pv);

    auto &layerParams =
        DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(
            layer_element);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams, "layer_material");
    }

    for (xml_coll_t ri(x_layer, _U(ring)); ri; ++ri) {
      xml_comp_t  x_ring   = ri;
      // double      r        = x_ring.r();
      //double      phi0     = x_ring.phi0(0);
      double      zstart   = x_ring.zstart();
      double      dz       = x_ring.dz(0);
      // int         nmodules = x_ring.nmodules();
      string      m_nam    = x_ring.moduleStr();
      Volume      m_vol    = modules[m_nam];
      //double      iphi     = 2 * M_PI / nmodules;
      //double      phi      = phi0;
      Placements& sensVols = sensitives[m_nam];

      // bool RSU_x_in_disk = true;
      // bool RSU_y_in_disk = true;
      // Starting values for first placement - these will change
      double RSU_x_pos = layer_x_offset;
      double RSU_y_pos = layer_rmin+(RSU_length/2);
      double RSU_z_pos = 0;
      //int RSU_id=0;
      string m_base;
      string zstring;
      double pos_multiplier[] = {-1,1};
      double tile_ymin = 1.5*RSU_width;
      double opening_ymin = 0;


      if (reflect){
	RSU_z_pos = -zstart - dz;
	zstring = "_neg";
      }
      else {
	RSU_z_pos = zstart + dz;
	zstring = "_pos";
      }

      // new module placement implementation
      //if (!reflect) {

	double layer_xmax = sqrt(layer_rmax*layer_rmax - (1.5*RSU_width + RSU_length)*(1.5*RSU_width + RSU_length));
	// Populate +ve x side of disk
	while (RSU_x_pos + 0.5*RSU_width < layer_xmax){
	  double layer_ymax = sqrt(layer_rmax*layer_rmax - (RSU_x_pos+0.5*RSU_width)*(RSU_x_pos+0.5*RSU_width));
	  while (RSU_y_pos + 0.5*RSU_length < layer_ymax){
	    for (auto multiplier : pos_multiplier){
	      m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_RSU%d");
	      DetElement module(layer_element, m_base + zstring, det_id);
	      pv = layer_vol.placeVolume(
					 m_vol, Transform3D(RotationZYX(0, 0, -M_PI / 2), Position(RSU_x_pos, multiplier*RSU_y_pos, RSU_z_pos)));
	      pv.addPhysVolID("module", mod_num);
	      module.setPlacement(pv);
	      for (size_t ic = 0; ic < sensVols.size(); ++ic) {
		PlacedVolume sens_pv = sensVols[ic];
		DetElement   comp_elt(module, sens_pv.volume().name(), mod_num);
		auto &comp_elt_params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
		comp_elt_params.set<string>("axis_definitions", "XZY");
		comp_elt.setPlacement(sens_pv);
		volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
	      }
	      ++mod_num;
	    }
	    RSU_y_pos += RSU_length;
	  }
	  RSU_x_pos += RSU_width;
	  RSU_y_pos = tile_ymin+0.5*RSU_length;
	  opening_ymin = 0;
	  if (RSU_x_pos - 0.5*RSU_width < layer_rmin + layer_x_offset){
	    opening_ymin = sqrt(layer_rmin*layer_rmin - (RSU_x_pos - 0.5*RSU_width - layer_x_offset)*(RSU_x_pos - 0.5*RSU_width - layer_x_offset));
	  }
	  if (opening_ymin > tile_ymin){
	    RSU_y_pos = opening_ymin+0.5*RSU_length;
	  }
	}
	// Reset positions to left of centre
	RSU_x_pos = layer_x_offset - RSU_width;
	RSU_y_pos = tile_ymin+0.5*RSU_length;
	RSU_y_pos = tile_ymin+0.5*RSU_length;
	opening_ymin = 0;
	if (RSU_x_pos + 0.5*RSU_width > -1*layer_rmin + layer_x_offset){
	  opening_ymin = sqrt(layer_rmin*layer_rmin - (RSU_x_pos + 0.5*RSU_width - layer_x_offset)*(RSU_x_pos + 0.5*RSU_width - layer_x_offset));
	}
	if (opening_ymin > tile_ymin){
	  RSU_y_pos = opening_ymin+0.5*RSU_length;
	}
	// Populate -ve x side of disk
	while (RSU_x_pos - 0.5*RSU_width > -1*layer_xmax){
	  double layer_ymax = sqrt(layer_rmax*layer_rmax - (RSU_x_pos-0.5*RSU_width)*(RSU_x_pos-0.5*RSU_width));
	  while (RSU_y_pos + 0.5*RSU_length < layer_ymax){
	    for (auto multiplier : pos_multiplier){
	      m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_RSU%d");
	      DetElement module(layer_element, m_base + zstring, det_id);
	      pv = layer_vol.placeVolume(
					 m_vol, Transform3D(RotationZYX(0, 0, -M_PI / 2), Position(RSU_x_pos, multiplier*RSU_y_pos, RSU_z_pos)));
	      pv.addPhysVolID("module", mod_num);
	      module.setPlacement(pv);
	      for (size_t ic = 0; ic < sensVols.size(); ++ic) {
		PlacedVolume sens_pv = sensVols[ic];
		DetElement   comp_elt(module, sens_pv.volume().name(), mod_num);
		auto &comp_elt_params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
		comp_elt_params.set<string>("axis_definitions", "XZY");
		comp_elt.setPlacement(sens_pv);
		volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
	      }
	      ++mod_num;
	    }
	    RSU_y_pos += RSU_length;
	  }
	  RSU_x_pos -= RSU_width;
	  RSU_y_pos = tile_ymin+0.5*RSU_length;
	  opening_ymin = 0;
	  if (RSU_x_pos + 0.5*RSU_width > -1*layer_rmin + layer_x_offset){
	    opening_ymin = sqrt(layer_rmin*layer_rmin - (RSU_x_pos + 0.5*RSU_width - layer_x_offset)*(RSU_x_pos + 0.5*RSU_width - layer_x_offset));
	  }
	  if (opening_ymin > tile_ymin){
	    RSU_y_pos = opening_ymin+0.5*RSU_length;
	  }
	}
	//layer_xmax = sqrt(layer_rmax*layer_rmax - (1.5*RSU_width)*(1.5*RSU_width));
	for (auto multiplier_x : pos_multiplier){
	  for (double multiplier_y : {-1, 0, 1}){
	    // Reset positions before placing horizontal bar
	    RSU_x_pos = 0;
	    RSU_y_pos = 0;
	    double n_rsu_x {0};
	    while (abs(RSU_x_pos) + 0.5*RSU_length < layer_rmax-RSU_length){
	      RSU_x_pos = layer_x_offset + multiplier_x*(sqrt(layer_rmin*layer_rmin - (multiplier_y*0.5*RSU_width)*multiplier_y*0.5*RSU_width) + (n_rsu_x+0.5)*RSU_length);
	      // RSU_x_pos = layer_x_offset + multiplier_x*(sqrt(layer_rmin*layer_rmin - (multiplier_y*0.5*RSU_width)*multiplier_y*0.5*RSU_width) + (0.5)*RSU_length);
	      RSU_y_pos = multiplier_y*RSU_width;
	      m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_RSU%d");
	      DetElement module(layer_element, m_base + zstring, det_id);
	      pv = layer_vol.placeVolume(
					 m_vol, Transform3D(RotationZYX(0, -M_PI / 2, -M_PI / 2), Position(RSU_x_pos, RSU_y_pos, RSU_z_pos)));
	      pv.addPhysVolID("module", mod_num);
	      module.setPlacement(pv);
	      for (size_t ic = 0; ic < sensVols.size(); ++ic) {
		PlacedVolume sens_pv = sensVols[ic];
		DetElement   comp_elt(module, sens_pv.volume().name(), mod_num);
		auto &comp_elt_params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
		comp_elt_params.set<string>("axis_definitions", "XZY");
		comp_elt.setPlacement(sens_pv);
		volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
	      }
	      ++mod_num;
	      ++n_rsu_x;
	    }
	  }
	}
    }
  }
  pv = motherVol.placeVolume(assembly, Position(0, 0, (reflect ? -1.0e-9 : 1.0e-9)));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(epic_TileEndcapTracker, create_detector)
