from dropsondes import dropsondes

ds = DropSondes()

def get_sounding(filepath):
	skewt_data = ds.get_sounding_data(filepath)
	
	return skewt_data
	
def get_dropsonde(filepath):
	skewt_data = ds.get_dropsonde_data(filepath)
	
	return skewt_data
	
		
	