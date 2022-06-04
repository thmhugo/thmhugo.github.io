function createLink(title, id) {
    var link = document.createElement("a");
    var linkText = document.createTextNode(title);

    link.appendChild(linkText);
    link.title = title;
    link.href = "#" + id;
    return link;
}

function referenceFigures() {
    var toRef = document.getElementsByClassName("autoref");
    var allCaptionedFigures = document.getElementsByClassName("caption");
    var figureNumbers = {};

    for (var i = 0; i < allCaptionedFigures.length; i++) {
        figureNumbers[allCaptionedFigures[i].id] = i + 1;
    }

    const l = toRef.length;
    for (var i = 0; i < l; i++) {
        var ref = toRef[0]; // Sucks but replaceWith consumes the array toRef.
        const figureId = ref.outerText;
        const title = "Figure " + figureNumbers[figureId];
        ref.replaceWith(createLink(title, figureId));
    }
}