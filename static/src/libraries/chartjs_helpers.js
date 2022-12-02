const updateData = (chart, data) => {
    chart.data.labels = data.map(entry => entry.mz);
    chart.data.datasets[0].data = data.map(entry => entry.abundance);
    chart.update();
}

export default updateData;