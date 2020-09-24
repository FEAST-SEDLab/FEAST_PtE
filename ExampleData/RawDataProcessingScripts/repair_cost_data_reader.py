def repair_cost_data_reader(file_in=None, notes=None, file_out=None):
    """
    reads a csv file containing repair cost data and stores the data as a standard python list. The first column of the
    input csv must contain the leak cost data. If the column has a header, it must be textual (not a number).
    file_path_in: Path to a data file
    notes:        Comments on the data
    file_out:     Path to store the data at
    """
    import csv
    from feast.input_data_classes import RepairData
    import pickle

    repair_costs = []

    with open(file_in) as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            if row[0][0].isdigit():
                repair_costs.append(float(row[0]))

    repair_out = RepairData(raw_file_name=file_in, notes=notes)
    repair_out.define_data(repair_costs=repair_costs)

    if file_out is not None:
        pickle.dump(repair_out, open(file_out, 'wb'))
    return repair_out
