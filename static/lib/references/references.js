function createLink(title, id) {
    var link = document.createElement("a");
    var linkText = document.createTextNode(title);

    link.appendChild(linkText);
    link.title = title;
    link.href = "#" + id;
    return link;
}

function indexFigures(figures) {
    figureNumbers = {};
    for (var i = 0; i < figures.length; i++) {
        figureNumbers[figures[i].id] = i + 1;
    }
    return figureNumbers;
}

function indexAlgorithms(algorithms) {
    var algorithmsLabels = {};

    for (var algo of algorithms)
    {
        if (algo.id !== "") // Just the labeled (with an id) ones
        {
            algorithmsLabels[algo.id] = algo.getElementsByClassName("ps-keyword")[0].innerHTML;
        }
    }

    return algorithmsLabels;
}

function makeReferences() {
    var toRef = document.getElementsByClassName("autoref");
    var allCaptionedFigures = document.getElementsByClassName("caption");
    var allCaptionedAlgorithms = document.getElementsByClassName("ps-root");


    var figureNumbers = indexFigures(allCaptionedFigures);
    var algorithmsLabels = indexAlgorithms(allCaptionedAlgorithms);

    const l = toRef.length;

    for (var i = 0; i < l; i++) {
        var ref = toRef[0]; // Sucks but replaceWith consumes the array toRef.

        const id = ref.outerText;
        var parent = document.getElementById(id).parentNode;

        if (parent.classList.contains("figure") || parent.classList.contains("figure-table-container"))
        {
            const title = "Figure " + figureNumbers[id];
            ref.replaceWith(createLink(title, id));
        }
        else
        {
            const title = algorithmsLabels[id];
            ref.replaceWith(createLink(title, id));
        }
    }
}